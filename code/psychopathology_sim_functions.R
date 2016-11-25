#item variance formula: 
#  var(y_i) = loading_i^2 * var(factor) + var(e_i).
#
#So, if we have loading^2 under a standardized factor model (var = 1.0), and
#  a target item variance, then we can compute the error variance as:
#  var(e_i) = var(y_i) - loadingi^2
#
#Example: factor loading of 0.9, but target item variance of 1.0
#  var(e_i) = 1.0 - .9^2 = .19
computeResidvar <- function(targetitemvar, floadings, fvar=1.0) { 
  targetitemvar - apply(floadings^2, 2, sum)*fvar 
}

#simple function to work with simsem to return the data with simulation results so that it can be extracted for graphs
savedat <- function(out, data) {
  data #just return the dataset
}

#function to add error correlation to model (theta) 
addErrorCor <- function(theta, vp, corval) {
  #theta is a named square residual covariance matrix. Minimally, diagonal must be populated with residual variances
  #varpair is a 2-element vector of the two variables to correlate
  #corval is the desired residual correlation between variables
  stopifnot(nrow(theta)==ncol(theta))
  stopifnot(rownames(theta)==colnames(theta))
  stopifnot(length(vp) == 2L)
  stopifnot(all(vp %in% rownames(theta)))
  
  theta[vp[1],vp[2]] <- theta[vp[2],vp[1]] <- corval*sqrt(theta[vp[1],vp[1]])*sqrt(theta[vp[2],vp[2]])
  return(theta)
}

#main worker function for generating lavaan syntax for CFA models
buildLavaanSyntax <- function(varnames, lambda, theta, psi, psistart=FALSE, thetastart=FALSE) {
  #lambda is an nfactors x nvariables loadings matrix (lambda)
  #theta is the error/residual covariance matrix for observed dependent variables
  #psi is the covariance matrix of latent variables (factors)
  stopifnot(inherits(lambda, "matrix"))
  stopifnot(length(varnames) == ncol(lambda))
  stopifnot(nrow(psi) == nrow(lambda)) #fvar should have the same number of elements as rows of loadings
  stopifnot(nrow(psi) == ncol(psi)) #psi should be square (covariance matrix)
  syntax=c()
  fitsyntax <- c() #syntax for fitting simulated replications (doesn't get anything for free)
  #setup factor structure #lambda is nfactors x nvariables
  for (i in 1:nrow(lambda)) {
    syntax <- c(syntax, paste0("f", i, " =~ ", paste(lambda[i, which(lambda[i,] != 0)], varnames[which(lambda[i,] != 0)], sep='*', collapse=" + ")))
    fitsyntax <- c(fitsyntax, paste0("f", i, " =~ NA*", paste(varnames[which(lambda[i,] != 0)], collapse=" + "))) #force first loading to be free by prepending NA (standardized solution)
  }
  
  #setup factor variance structure (psi)
  for (i in 1:nrow(psi)) {
    for (j in 1:ncol(psi)) {
      if (i > j) next #only lower triangle 
      if (psi[i,j] != 0) { 
        #only add syntax if non-zero
        syntax <- c(syntax, paste0("f", i, " ~~ ", psi[i,j], "*", "f", j))
        if (psistart) { fitsyntax <- c(fitsyntax, paste0("f", i, " ~~ start(", psi[i,j], ")*", "f", j)) } #give it the right starting value
      }
      if (psi[i,j] == 0) { fitsyntax <- c(fitsyntax, paste0("f", i, " ~~ ", psi[i,j], "*", "f", j)) } #add explicit zero covariance to psi in fitting (identification for bifactor)       
      if (i==j) { fitsyntax <- c(fitsyntax, paste0("f", i, " ~~ ", psi[i,j], "*", "f", j)) } #add standardized factor variances (fixed)
    }
  }
  
  #setup error variance structure (theta)
  for (i in 1:nrow(theta)) {
    for (j in 1:ncol(theta)) {
      if (i > j) next #only lower triangle 
      if (theta[i,j] != 0) { #only add syntax if non-zero 
        syntax <- c(syntax, paste0(varnames[i], " ~~ ", theta[i,j], "*", varnames[j]))
        if (thetastart) {
          #adding variances to fit syntax for the moment (testing)
          fitsyntax <- c(fitsyntax, paste0(varnames[i], " ~~ start(", theta[i,j], ")*", varnames[j]))
        } else {
          #adding variances to fit syntax for the moment (testing)
          fitsyntax <- c(fitsyntax, paste0(varnames[i], " ~~ ", varnames[j]))
        }
      }
    }
  }
  return(list(simsyntax=paste(syntax, collapse="\n"), fitsyntax=paste(fitsyntax, collapse="\n")))
}

simCFAGraphs <- function(model, nreplications, n, graphmethods=c("EBICglasso", "pcor"), ...) { #, "IsingFit"
  require(simsem)
  require(parallel)
  require(abind)
  require(bootnet)
  syntax <- buildLavaanSyntax(model$varnames, model$lambda, model$theta, model$psi, ...) #build lavaan syntax for simulation

  #note that by default, the fitted models are provided only with a configural model on which free parameters are estimated
  simStruct <- simsem::sim(nRep=nreplications, model=syntax$fitsyntax, n=n, generate=syntax$simsyntax, 
      lavaanfun = "cfa", outfundata=savedat, multicore=FALSE) #TRUE)
  
  #obtain raw simulated data for network analysis
  dlist <- getExtraOutput(simStruct) #only seems to work with convergence of estimated model?
  #dlist <- simStruct@extraOut #manual (always works)
  
  #setup output structure
  outstruct <- list()
  outstruct$specification <- list(model=model, nreplications=nreplications, n=n, graphmethods=graphmethods, syntax=syntax) #snapshot of simulation structure
  outstruct$simsemout <- simStruct
  outstruct$simdata <- dlist
  
  #basic adjacency matrix based on pairwise pearson correlation (not partial)
  outstruct$adjmats[["pearson"]]$concat <- do.call(abind, list(lapply(outstruct$simdata, function(df) { cor(df) }), along=0))
  outstruct$adjmats[["pearson"]]$average <- apply(outstruct$adjmats[["pearson"]]$concat, c(2,3), mean)
  
  cl <- makeCluster(4) #defaults to PSOCK cluster
  clusterExport(cl, c("estimateNetwork"))#, "cd4.mle"))
  clusterEvalQ(cl, library(bootnet))
  clusterEvalQ(cl, library(qgraph))
  
  on.exit(try(stopCluster(cl)))
  
  for (method in graphmethods) {
    glist <- parLapply(cl, dlist, function(df) { estimateNetwork(df, default=method) })
    
    outstruct$graphs[[method]] <- glist #save raw graphs into output
    
    #concatenated and average adjacency matrices
    outstruct$adjmats[[method]]$concat <- do.call(abind, list(lapply(glist, "[[", "graph"), along=0))
    outstruct$adjmats[[method]]$average <- apply(outstruct$adjmats[[method]]$concat, c(2,3), mean)
    outstruct$adjmats[[method]]$se <- apply(outstruct$adjmats[[method]]$concat, c(2,3), plotrix::std.error)
    
    #compute centrality metrics (just the handful supported in bootnet)
    Long <- parLapply(cl, 1:length(glist), function(g) { 
          df <- centralityTable(glist[[g]], standardized=TRUE)
          df$graphNum <- g; return(df) })
    LongAll <- do.call(rbind, Long) #includes metrics for all replications (for graphing means and variation/uncertainty
    
    #browser()
    
    #manual calculation using igraph
    #ct = centralityTable(glist[[1]], standardized=FALSE)
    
    #gg <- igraph::graph_from_adjacency_matrix(outstruct$adjmats[["EBICglasso"]]$concat[1,,], mode="undirected", weighted=TRUE, diag=FALSE)
    #igraph::closeness(gg, normalized=TRUE) #blowing up igraph for unknown reasons
    
    #qgraph uses 1/correlation as weight, whereas igraph has to be specifically asked for this 
    #igraph::betweenness(gg, normalized=FALSE, weights=1/igraph::E(gg)$weight) #, weights=NULL)
    #centrality_auto(glist[[1]])$node.centrality
    
    #igraph::strength(gg)
    #igraph::closeness(gg, normalized=FALSE, weights=1/igraph::E(gg)$weight)
    #igraph::betweenness(gg, normalized=FALSE, weights=igraph::E(gg)$weight)
    #cor(igraph::closeness(gg, normalized=FALSE), unlist(subset(ct, measure=="Closeness", select=value)))
    #gg <- graph_from_adjacency_matrix(outstruct$adjmats[["pearson"]]$concat[1,,], mode="undirected", weighted=TRUE, diag=FALSE)
    
    #igraph::closeness(gg)
    
    #for some reason, pcor drops the "y" prefix
    #this will only work if the variables are named y<...>
    if (method == "pcor") {
      message("Adding missing 'y' prefix to node names for pcor")
      LongAll$node <- paste0("y", LongAll$node) 
    }
    
    browser()
    loadmaster <- data.frame(poploading=gdata::unmatrix(model$lambda, byrow=TRUE), 
        factor=paste0(rep(dimnames(model$psi)[[1]], each=ncol(model$lambda)), "_poploading"), node=model$varnames, stringsAsFactors=FALSE)
    
    #convert loadings into wide format (columns named by factors)
    loadmaster <- loadmaster %>% tidyr::spread(key=factor, value=poploading)
    
    LongAll <- inner_join(LongAll, loadmaster, by="node")
        
    #fitted_loadings <- select(sims[["simsemout"]]@coef, starts_with("f1=~"))
    fitted_loadings <- select(outstruct[["simsemout"]]@coef, matches("f\\d+=~"))
    fitted_loadings$graphNum <- 1:nrow(fitted_loadings)
    #names(fitted_loadings) <- sub("f1=~", "", names(fitted_loadings), fixed=TRUE)
    wide_fitted <- tidyr::gather(fitted_loadings, key="node", value="fittedloading", -graphNum) %>%
        tidyr::separate(node, into=c("factor", "node"), sep="=~") %>% mutate(factor=paste0(factor, "_fittedloading")) %>%
        tidyr::spread(key=factor, value=fittedloading)
    
    LongAll <- inner_join(LongAll, wide_fitted, by=c("node", "graphNum"))
    outstruct$gmetrics[[method]] <- LongAll #save centrality metrics into output
    
    outstruct[["graph_v_factor"]][[method]]$metric_v_loadings <- LongAll %>% 
        gather(key=factor, value=loading, matches("f\\d+.*loading")) %>% # f1_poploading, f2_poploading, f1_fittedloading, f2_fittedloading) %>%
        separate(col=factor, into=c("factor", "loadingtype")) %>% 
        spread(key="loadingtype", value="loading") %>% filter(!is.na(poploading) & poploading != 0)

    gg <- igraph::graph_from_adjacency_matrix(glist[[1]]$graph, mode="undirected", weighted=TRUE, diag=FALSE)
    gg <- delete.edges(gg, which(E(gg)$weight < 0))
    
    ct <- centralityTable(glist[[1]], standardized=FALSE)
    #ct <- centralityTable(gg, standardized=FALSE)
    cor(igraph::strength(gg), subset(ct, measure=="Strength", select=value))
    cor(igraph::closeness(gg, weights=1/E(gg)$weight), subset(ct, measure=="Closeness", select=value), use="pairwise.complete.obs")
    cor(igraph::betweenness(gg, weights=1/E(gg)$weight), subset(ct, measure=="Betweenness", select=value), use="pairwise.complete.obs")
    centrality_auto(glist[[1]])
    centrality_auto(glist[[1]])
    
    
    adjmat <- glist[[1]]$graph
    weights <- adjmat[lower.tri(adjmat)]
    length(weights[weights != 0])
    igraph::closeness(gg, normalized=TRUE, weights=1/igraph::E(gg)$weight) #blowing up igraph for unknown reasons
    igraph::betweenness(gg, normalized=FALSE, weights=1/E(gg)$weight) 
    #qgraph uses 1/correlation as weight, whereas igraph has to be specifically asked for this 
    #igraph::betweenness(gg, normalized=FALSE, weights=1/igraph::E(gg)$weight) #, weights=NULL)
    #centrality_auto(glist[[1]])$node.centrality
    
    cmat <- cor_auto(dlist[[3]], detectOrdinal = FALSE)
    EBICgraph <- qgraph(cmat, graph = "glasso", sampleSize = nrow(dlist[[1]]),
        tuning = 0.5, layout = "spring", title = "BIC", details = TRUE)
    
    #this is the estimation function under the hood of estimateNetwork
    tt = EBICglasso(S=cmat, n=nrow(dlist[[1]]), gamma=0.5, penalize.diagonal=FALSE, nlambda=100,
        lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE, 
        countDiagonal = FALSE, refit = FALSE) #penalizeMatrix, 
    
    library(qgraph)
    fromq <- as.igraph(EBICgraph)
    s1 <- sort(E(fromq)$weight)
    gg <- igraph::graph_from_adjacency_matrix(tt, mode="undirected", weighted=TRUE, diag=FALSE)
    s2 <- sort(E(gg)$weight)
    all.equal(s1, s2) #yes, these match.
    
    #for sanity, let's switch to igraph for everything after the network is estimated...
    #gg <- delete.edges(gg, which(E(gg)$weight < 0))
    E(gg)$weight <- abs(E(gg)$weight)
    closeness(gg, weights=1/E(gg)$weight)
    cor(betweenness(gg, weights=1/E(gg)$weight), centrality(EBICgraph, pkg="igraph")$Betweenness) #looks like qgraph is treating as directed? all values doubled wrt igraph

    #so, this matches. CentralityTable standardizes (z-scores) all measures... unclear why closeness blows up
    #look at the examples that are bad. dlist[[3]] in this iteration
    cor(closeness(gg, weights=1/E(gg)$weight), centrality(EBICgraph, pkg="igraph")$Closeness)
    betweenness(gg)
    
    
    g2 <- estimateNetwork(dlist[[3]], default="EBICglasso")
    identical(g2$graph, tt) #TRUE
    
    
    outstruct$graph_v_factor[[method]]$corr_v_fitted <- outstruct$graph_v_factor[[method]]$metric_v_loadings %>%
        group_by(graphNum, measure, factor) %>% do({
              if (all(is.na(.$value))) { browser("missing graph") }
              tryCatch(sd(.$fittedloading) < .01 || sd(.$value), error=function(x) { browser()})
              if (sd(.$fittedloading) < .01 || sd(.$value) < .01) { data.frame(cv=NA)
              } else { data.frame(graphNum=.$graphNum[1], measure=.$measure[1], cv=cor(.$fittedloading, .$value)) } }) %>% #summarize(cv=cor(loading,value)) %>%
        group_by(measure, factor) %>% summarize(mcv=mean(cv, na.rm=TRUE))
    
  }
  
  
  
  return(outstruct)
}

##First demonstration: correlation of one-factor CFA loadings with centrality measures
##a function to simulate multiple one-factor CFAs:
## 10, 15, or 20 indicators
## a basis for sampling loadings like uniform from seq(.4, .95, .05)
## a number of replications for each setting.

demo1 <- function(nexamples=2, nindicators=10, nreplications=200, n=400, loadingsampler=seq(.4, .95, .05), 
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL) {
  # nexamples is the number of population factor models from which data are simulated and fit.
  #   Each example draws a set of population factor loadings from the loadingsampler.
  # nindicators is the number of indicators of the factor
  # nreplications is the number of datasets drawn from the model
  # n is the sample size
  # loadingsampler is the set of population factor loadings that are drawn uniformly with replacement
  # errorvars specifies whether to have equal item residual variance (eqresidvar) or equal observed variance
  # ivar is the itevm variance (or residual variance)
  #
  require(dplyr)
  require(tidyr)
  require(matrixcalc)
  
  lapply(1:nexamples, function(i) {
        lambda <- as.matrix(rbind(sample(loadingsampler, nindicators, replace=TRUE)))
        varnames <- paste0("y", 1:ncol(lambda))
        fvar <- 1.0 #factor variance (standardized)
        if (errorvars=="eqresidvar") {
          errorvars <- rep(ivar, length(varnames)) #fixed error specification          
        } else if (errorvars=="eqvars") {
          errorvars <- computeResidvar(targetitemvar=ivar, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances          
        }
        
        theta <- diag(as.vector(errorvars)) #0 resid cov by default
        rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work
        
        #thetacorlist contains an edge list specifying the residual correlation of two items
        if (!is.null(thetacorlist)) {
          for (pair in thetacorlist) {
            theta <- addErrorCor(theta, unlist(pair[c(1,2)]), pair[[3]]) #positions 1 and 2 are the named nodes/items, 3 is the target correlation  
          }
        }
        
        if (!is.positive.definite(theta)) {
          message("theta (error cov matrix) is not positive definite. Cannot continue")
          print(theta)
          return(NULL)
        }
        
        psi <- diag(fvar) #zero covariance in factor structure at the moment
        dimnames(psi) <- list(f=paste0("f", 1:nrow(psi)), f=paste0("f", 1:ncol(psi)))
        
        #model specification structure. Currently just wrapping up individual arguments above into a list
        m <- list(
            varnames=varnames, #vector of variable names
            lambda=lambda, #nitems x nfactors loadings matrix
            theta=theta, #covariance (residual) matrix for observed items
            psi=psi #covariance matrix for factors
        )
        
        sims <- simCFAGraphs(m, nreplications=nreplications, n=n)
                
        return(sims)
      })
}


demomaster <- function(nexamples=2, nindicators=20, nfactors=2, nreplications=200, n=400, 
    loadings="random", loadingsampler=seq(.4, .95, .05),
    prop_residcor=0.0, mag_residcor=function() { runif(1, 0.3, 0.7) },
    prop_crossload=0.0, mag_crossload=function() { runif(1, 0.3, 0.7) },
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL) {
  
  # - nexamples: the number of population factor models from which data are simulated and fit.
  #    Each example draws a set of population factor loadings from the loadingsampler.
  # - nindicators: the number of total indicators in the model
  # - nfactors: the number of factors in the model. By default, indicators will be distributed equally across factors
  # - nreplications: the number of datasets drawn from the model-implied mean and covariance structure
  # - n: the sample size of each replication
  # - loadings: the factor loadings approach for simulation. The following conditions are supported:
  #      1) A scalar value (e.g., 0.8). In this case, this value will be the loading for all indicators on their respective factor
  #      2) "random". This specifies to draw factor loadings randomly from loadingsampler for each example (but fixed across replications within example)
  #      3) An nfactors x nindicators matrix. This is a custom lambda (factor loadings) matrix and overrides nfactors and nindicators.
  # - loadingsampler: the set of population factor loadings that are drawn uniformly with replacement for each example when loadings="random"
  # - prop_residcor: the proportion of indicators that have a random residual correlation
  # - mag_residcor: a function that returns a single random number for a residual correlation
  # - prop_crossload: the proportion of indicators that are allowed to cross-load randomly on a secondary factor
  # - mag_crossload: a function that returns a single random number for a secondary factor loading
  # - errorvars specifies whether to have equal item residual variance (eqresidvar) or equal observed variance
  # - ivar is the item variance (or residual variance)
  # - thetacorlist is an edge-list style custom error cor specification. Each element is a three-element list of the form: list("y1", "y2", 0.5). 

  require(dplyr)
  require(tidyr)
  require(matrixcalc)

  randomloadings <- TRUE
  
  if (is.matrix(loadings)) {
    message("Using custom loadings matrix (lambda) for all simulations")
    lambda <- loadings
    nfactors <- nrow(lambda)
    nindicators <- ncol(lambda)
    randomloadings <- FALSE
  } else {
    #number of indicators must be a multiple of number of factors (equal number of indicators per factor)
    stopifnot (nindicators %% nfactors == 0)
    perfactor <- nindicators / nfactors

    lambda <- matrix(0, nrow=nfactors, ncol=nindicators)
    
    if (is.numeric(loadings[1L]) && length(loadings) == 1L) {
      message("Using a fixed factor loading for all factors: ", loadings)
      randomloadings <- FALSE
      for (i in 1:nfactors) {
        lambda[i,(i-1)*perfactor + 1:perfactor] <- loadings
      }
    }
  }
  
  varnames <- paste0("y", 1:ncol(lambda))
  dimnames(lambda) <- list(factor=paste0("f", 1:nrow(lambda)), indicator=varnames)
  
  theta <- matrix(0, nrow=nindicators, ncol=nindicators)
  rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work
  
  if (prop_residcor > 0) {
    lt <- which(lower.tri(theta), arr.ind=TRUE) #positions of lower diagonal in rows, columns
    lt_names <- apply(lt, c(1,2), function(x) { varnames[x] }) #convert to varnames (addErrorCorr requests this)
    rc_positions <- lt_names[sample(1:nrow(lt_names), floor(prop_residcor*nindicators), replace=FALSE),] #consistent positions of error correlations across examples (easier to track)    
  } else { rc_positions <- matrix(numeric(0), nrow=0, ncol=2) }
  
  if (prop_crossload > 0 && nfactors > 1) {
    fpos <- as.data.frame(which(lambda > 0, arr.ind=TRUE), row.names=1:ncol(lambda)) #positions of non-zero factor loadings (need to be careful of rownames failing unless specified)
    
    #rule out any indicators that already have cross-loadings from the random cross-loading process
    loadtab <- as.data.frame(table(fpos[,2]))
    nocross <- with(loadtab, as.numeric(Var1[Freq == 1]))
    fpos <- subset(fpos, indicator %in% nocross)

    fnums <- unique(fpos$factor)
    cl_positions <- fpos[sample(1:nrow(fpos), floor(prop_crossload*nindicators), replace=FALSE),] #consistent positions of cross-loadings across examples (easier to track)

    #cl_positions now contains random samples of indicators without cross-loadings. But we need to actually cross-load it onto an eligible factor
    for (i in 1:nrow(cl_positions)) {
      elig <- fnums[fnums != cl_positions[i,"factor"]]
      cl_positions[i,"factor"] <- elig[sample(length(elig), 1)] #http://stackoverflow.com/questions/13990125/sampling-in-r-from-vector-of-varying-length 
    }
  } else { cl_positions <- matrix(numeric(0), nrow=0, ncol=2) }

  #generate examples
  lapply(1:nexamples, function(i) {
        if (randomloadings) {
          for (i in 1:nfactors) {
            lambda[i,(i-1)*perfactor + 1:perfactor] <- sample(loadingsampler, perfactor, replace=TRUE)
          }      
        }

        fvar <- rep(1.0, nfactors) #standardized factor variance
        psi <- diag(fvar) #zero covariance in factor structure at the moment
        dimnames(psi) <- list(f=paste0("f", 1:nrow(psi)), f=paste0("f", 1:ncol(psi)))
        
        #allow for random cross-loadings
        if (nrow(cl_positions) > 0) {
          for (i in 1:nrow(cl_positions)) {
            lambda[cl_positions[i,1], cl_positions[i,2]] <- mag_crossload()
          }
        }
      
        if (errorvars=="eqresidvar") {
          errorvars <- rep(ivar, length(varnames)) #fixed error specification          
        } else if (errorvars=="eqvars") {
          errorvars <- computeResidvar(targetitemvar=ivar, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances          
        }
        
        diag(theta) <- as.vector(errorvars) #0 resid cov by default

        #loop over random residual correlations to incorporate
        if (prop_residcor > 0) {
          for (i in 1:nrow(rc_positions)) {
            theta <- addErrorCor(theta, rc_positions[i,], mag_residcor()) #positions 1 and 2 are the numeric nodes/items, 3 is the target correlation  
          }
        } 
                
        #thetacorlist contains an edge list specifying the residual correlation of two items
        #note that this could override the random residual correlations from the previous step
        if (!is.null(thetacorlist)) {
          for (pair in thetacorlist) {
            theta <- addErrorCor(theta, unlist(pair[c(1,2)]), pair[[3]]) #positions 1 and 2 are the named nodes/items, 3 is the target correlation  
          }
        }
        
        #check that theta has not become degenerate
        if (!is.positive.definite(theta)) {
          message("theta (error cov matrix) is not positive definite. Cannot continue")
          print(round(theta, 2))
          return(NULL)
        }
        
        #model specification structure. Currently just wrapping up individual arguments above into a list
        m <- list(
            varnames=varnames, #vector of variable names
            lambda=lambda, #nitems x nfactors loadings matrix
            theta=theta, #covariance (residual) matrix for observed items
            psi=psi #covariance matrix for factors
        )
        
        message("About to run the following model: ")
        print(m)
        
        sims <- simCFAGraphs(m, nreplications=nreplications, n=n)
        
        return(sims)
      })
}