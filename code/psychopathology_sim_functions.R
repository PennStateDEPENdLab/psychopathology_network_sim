#item variance formula: 
#  var(y_i) = loading_i^2 * var(factor) + var(e_i).
#
#So, if we have loading^2 under a standardized factor model (var = 1.0), and
#  a target item variance, then we can compute the error variance as:
#  var(e_i) = var(y_i) - loadingi^2
#
#Example: factor loading of 0.9, but target item variance of 1.0
#  var(e_i) = 1.0 - .9^2 = .19
#
#Also note that item reliability is equal to the squared loading (explained variance)
#  when item variance = 1.0 and factor variance = 1.0 (see Muthen & Muthen 2002 SEM Power study)
#  Item reliability = (lambda^2*psi)/(lambda^2*psi + theta)
#  where lambda is the item's factor loading, psi is factor variance, and theta is item residual variance
computeResidvar <- function(targetitemvar, floadings, fvar=1.0) {
  #fvar is a vector with one variance per factor
  #need to apply appropriate variance to its indicators
  #use the sweep function to row-wise multiply by the correct element of fvar
  targetitemvar - apply(sweep(floadings^2, 1, fvar, "*"), 2, sum)
  
  #targetitemvar - apply(floadings^2, 2, sum)*fvar #this will not work as expected for multi-factor models  
}

#simple function to work with simsem to return the data with simulation results so that it can be extracted for graphs
savedat <- function(out, data) {
  list(lvobj=out, data=data) #return the dataset and lavaan object
  #data #return the dataset
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
buildLavaanSyntax <- function(varnames, lambda, lambdaconstraint=NULL, theta, psi, psistart=FALSE, thetastart=FALSE) {
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
    nzlambda <- which(lambda[i,] != 0)
    
    if (!is.null(lambdaconstraint)) {
      thiscon <- sapply(lambdaconstraint[i,], function(x) { if(x==0) "NA" else paste0("v", x) })
      whichcon <- which(thiscon != "NA") #position of variables having a constraint
      
      #trap condition where first loading on factor is part of an equality constraint
      #in this case, we need to specify that indicator twice, one for freeing loading NA*y1 and once for constraint v1*y1
      #see modifiers section here: http://lavaan.ugent.be/tutorial/syntax2.html
      if (thiscon[nzlambda][1] != "NA") {
        prepend <- paste0("NA*", varnames[nzlambda][1], " + ")
      } else { prepend <- "" }
      
      if (length(whichcon) > 0L) { append <- paste0(" + ", paste(thiscon[whichcon], varnames[whichcon], sep='*', collapse=" + ")) } else { append <- "" }
      syntax <- c(syntax, paste0("f", i, " =~ ", paste(lambda[i, nzlambda], varnames[nzlambda], sep='*', collapse=" + "), append)) #append constraints

      fitsyntax <- c(fitsyntax, paste0("f", i, " =~ ", prepend, paste(thiscon[nzlambda], varnames[nzlambda], sep="*", collapse=" + "))) #force first loading to be free by prepending NA (standardized solution)  
    } else {
      syntax <- c(syntax, paste0("f", i, " =~ ", paste(lambda[i, nzlambda], varnames[nzlambda], sep='*', collapse=" + ")))
      fitsyntax <- c(fitsyntax, paste0("f", i, " =~ NA*", paste(varnames[nzlambda], collapse=" + "))) #force first loading to be free by prepending NA (standardized solution)  
    }
    
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

simCFAGraphs <- function(model, nreplications, n, parallel=4, graphmethods=c("EBICglasso", "pcor", "cor.shrink"), saveLavObj=FALSE, saveSimsemout=FALSE, seed=NULL, ...) { #, "IsingFit"
  require(simsem)
  require(parallel)
  require(abind)
  require(bootnet)
  require(corpcor)
  require(doParallel)
  require(igraph)
  require(tidyr)
  syntax <- buildLavaanSyntax(varnames=model$varnames, lambda=model$lambda, lambdaconstraint=model$lambdaconstraint, 
                              theta=model$theta, psi=model$psi, ...) #build lavaan syntax for simulation
  #note that by default, the fitted models are provided only with a configural model on which free parameters are estimated
  
  #randomize seed
  if (is.null(seed)) {
    seed = as.POSIXlt(Sys.time())
    seed = 1000*(seed$hour*3600 + seed$min*60 + seed$sec)
  }

  simStruct <- simsem::sim(nRep=nreplications, model=syntax$fitsyntax, n=n, generate=syntax$simsyntax, 
      lavaanfun = "cfa", outfundata=savedat, multicore=FALSE, seed=seed) #TRUE)
  
  #obtain raw simulated data for network analysis
  dlist <- getExtraOutput(simStruct) #only seems to work with convergence of estimated model?
  #dlist <- simStruct@extraOut #manual (always works)
  
  #setup output structure
  outstruct <- list()
  outstruct$specification <- list(model=model, nreplications=nreplications, n=n, graphmethods=graphmethods, syntax=syntax) #snapshot of simulation structure
  outstruct$simsemout <- simStruct
  outstruct$simdata <- lapply(dlist, "[[", "data")
  if (saveLavObj) {
    outstruct$lvobj <- lapply(dlist, "[[", "lvobj")
  }
  
  #basic adjacency matrix based on pairwise pearson correlation (not partial)
  outstruct$adjmats[["pearson"]]$concat <- do.call(abind, list(lapply(outstruct$simdata, function(df) { cor(df) }), along=0))
  outstruct$adjmats[["pearson"]]$average <- apply(outstruct$adjmats[["pearson"]]$concat, c(2,3), mean)
  
  #switch away from parLapply because it doesn't easily allow for serial execution, which is important if parallelism is in outer loop
  if (!is.null(parallel) && parallel > 0) {
    cl <- makeCluster(parallel) #defaults to PSOCK cluster
    #clusterExport(cl, c("estimateNetwork"))#, "cd4.mle"))
    #these were necessary under parLapply, but not under the new foreach using igraph directly
    #clusterExport(cl, c("cor_auto", "cor2pcor", "EBICglasso")) #should be handled by foreach
    #clusterEvalQ(cl, library(bootnet))
    #clusterEvalQ(cl, library(qgraph))
    #clusterEvalQ(cl, library(igraph))
    #clusterEvalQ(cl, library(tidyr))
    
    registerDoParallel(cl)
    
    on.exit(try(stopCluster(cl)))
    `%op%` <- `%dopar%` #https://stat.ethz.ch/pipermail/r-sig-hpc/2013-January/001575.html
  } else {
    #removed the doseq backend so that a parent function can be made parallel without overriding the registered backend
    #registerDoSEQ() #formally register sequential backend
    `%op%` <- `%do%`
  }
  
  for (method in graphmethods) {
    #achieve good speedup here using parallel execution
    glist <- foreach(df=outstruct$simdata, .packages=c("qgraph", "igraph", "tidyr", "corpcor")) %op% {
          if (method=="EBICglasso") {
            cmat <- cor_auto(df, detectOrdinal = TRUE, ordinalLevelMax = 7, missing = "pairwise")
            adjmat <- EBICglasso(S=cmat, n=nrow(df), gamma=0.5, penalize.diagonal=FALSE, nlambda=100,
                lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = TRUE, 
                countDiagonal = FALSE, refit = FALSE)
          } else if (method=="pcor") {
            cmat <- cor_auto(df, detectOrdinal = TRUE, ordinalLevelMax = 7, missing = "pairwise")
            adjmat <- cor2pcor(cmat)
            dimnames(adjmat) <- list(names(df), names(df)) #function drops the names
          } else if (method=="cor.shrink") { #shrinkage estimator of marginal association
            adjmat <- corpcor::cor.shrink(df)
            dimnames(adjmat) <- list(names(df), names(df)) #function drops the names
          }
          gg <- igraph::graph_from_adjacency_matrix(adjmat, mode="undirected", weighted=TRUE, diag=FALSE)
          return(gg)
          
          #estimateNetwork(df, default=method)
        }

    outstruct$graphs[[method]] <- glist #save raw graphs into output

    #concatenated and average adjacency matrices
    #use g[] notation to convert igraph to weighted matrix
    outstruct$adjmats[[method]]$concat <- do.call(abind, list(lapply(glist, function(g) { as.matrix(g[]) }), along=0))
    outstruct$adjmats[[method]]$average <- apply(outstruct$adjmats[[method]]$concat, c(2,3), mean)
    outstruct$adjmats[[method]]$se <- apply(outstruct$adjmats[[method]]$concat, c(2,3), plotrix::std.error)
    
    #compute centrality measures for each graph (just closeness, betweenness, and strength for now to match qgraph)
    #parLapply does not help with speedup here (probably because igraph code is optimized and compiled)
    Long <- lapply(1:length(glist), function(i) {
          gg <- glist[[i]] #igraph object
          E(gg)$weight <- abs(E(gg)$weight) #this is how qgraph solves negative edges (adopt for now)

          ct <- data.frame(graphNum=i, node=V(gg)$name, 
              closeness=closeness(gg, weights=1/E(gg)$weight, normalized=TRUE), #1/weight transformation because based on path-based distance metric
              betweenness=betweenness(gg, weights=1/E(gg)$weight), #1/weight transformation to treat as distance function
              strength=strength(gg), stringsAsFactors=FALSE)
          ct %>% gather(measure, value, closeness, betweenness, strength)
        })
    
    #compute centrality metrics (just the handful supported in bootnet)
    #Long <- parLapply(cl, 1:length(glist), function(g) { 
    #      df <- centralityTable(glist[[g]], standardized=FALSE)
    #      df$graphNum <- g; return(df) })

    LongAll <- do.call(rbind, Long) #includes metrics for all replications (for graphing means and variation/uncertainty
    
    #cor2pcor drops dimnames of matrix. For estimateNetwork, need to add back "y" prefix
    #this will only work if the variables are named y<...>
    #if (method == "pcor") {
    #  message("Adding missing 'y' prefix to node names for pcor")
    #  LongAll$node <- paste0("y", LongAll$node) 
    #}
    
    loadmaster <- data.frame(poploading=gdata::unmatrix(model$lambda, byrow=TRUE), 
        factor=paste0(rep(dimnames(model$psi)[[1]], each=ncol(model$lambda)), "_poploading"), 
        node=model$varnames, stringsAsFactors=FALSE)
    
    #convert loadings into wide format (columns named by factors)
    loadmaster <- loadmaster %>% tidyr::spread(key=factor, value=poploading)

    LongAll <- inner_join(LongAll, loadmaster, by="node")

    #need to trim out equality constraint baggage from coefficient names
    #example: v1 <- (f3=~y2)
    names(outstruct[["simsemout"]]@coef) <- sub("\\w+ <- \\(([^)]+)\\)", "\\1", names(outstruct[["simsemout"]]@coef), perl=TRUE)

    fitted_loadings <- select(outstruct[["simsemout"]]@coef, matches("f\\d+=~"))
    fitted_loadings$graphNum <- 1:nrow(fitted_loadings)

    wide_fitted <- tidyr::gather(fitted_loadings, key="node", value="fittedloading", -graphNum) %>%
        tidyr::separate(node, into=c("factor", "node"), sep="=~") %>% mutate(factor=paste0(factor, "_fittedloading")) %>%
        tidyr::spread(key=factor, value=fittedloading)
    
    LongAll <- inner_join(LongAll, wide_fitted, by=c("node", "graphNum"))
    
    #merge population and fitted residual correlations in theta to the centrality data.frame
    #crazy regexp is to remove any constraint coefficients, labeled [ <conname> ] - [ <conname> ]
    coefnames <- strsplit(grep("\\s*\\[\\s*\\w+\\s*\\]\\s*-\\s*\\[\\s*\\w+\\s*\\]\\s*", names(outstruct[["simsemout"]]@coef), perl=TRUE, value=TRUE, invert=TRUE), "~~|=~")
    
    stopifnot(all(sapply(coefnames, length) == 2))
    coefnames <- data.frame(do.call(rbind, coefnames), stringsAsFactors=FALSE)
    rpos <- with(coefnames, which(substr(X1, 1, 1)=="y" & substr(X2, 1, 1)=="y" & X1 != X2)) #filter(residcov, substr(X1, 1, 1)=="y" & substr(X2, 1, 1)=="y" & X1 != X2) 

    if (length(rpos) > 0L) {
      #populate residual covariance structure
      residcov <- outstruct[["simsemout"]]@coef[,rpos, drop=F]
      residlist <- list()
      
      #want to have a structure that looks like
      #graphNum node rvar rcorr_pop rcorr_fitted
      #       1   y1   y9      0.55         0.52
      #       2   y1   y9      0.55         0.59
      #       1   y9   y1      0.55         0.52
      #       2   y9   y1      0.55         0.59
      
      for (i in 1:length(rpos)) {
        pv <- model$theta[coefnames[rpos[i],1], coefnames[rpos[i],2]] #population value
        residlist[[i]] <- rbind(
            data.frame(node=coefnames[rpos[i],1], rvar=coefnames[rpos[i],2], rcorr_pop=pv, rcorr_fitted=residcov[,i], graphNum=1:nrow(residcov), stringsAsFactors=FALSE),
            data.frame(node=coefnames[rpos[i],2], rvar=coefnames[rpos[i],1], rcorr_pop=pv, rcorr_fitted=residcov[,i], graphNum=1:nrow(residcov), stringsAsFactors=FALSE)
        )
      }
      residcov <- do.call(rbind, residlist)
      
      #left join because residcov will only affect a select few nodes, NAs populated elsewhere
      LongAll <- left_join(LongAll, residcov, by=c("node", "graphNum"))  
    }
    
    outstruct$gmetrics[[method]] <- LongAll #save centrality metrics into output with population and fitted loadings
    
    #align centrality measures with factor loadings
    outstruct[["graph_v_factor"]][[method]]$metric_v_loadings <- LongAll %>% 
        gather(key=factor, value=loading, matches("f\\d+.*loading")) %>% # f1_poploading, f2_poploading, f1_fittedloading, f2_fittedloading) %>%
        separate(col=factor, into=c("factor", "loadingtype")) %>% 
        spread(key="loadingtype", value="loading") %>% filter(!is.na(poploading) & poploading != 0)

    #debugging
    #outstruct$graph_v_factor[[method]]$metric_v_loadings %>%
    #    group_by(graphNum, measure, factor) %>% summarize(centsd=sd(value), loadsd=sd(fittedloading)) %>% 
    #    filter(measure=="closeness") %>% arrange(desc(centsd))
    
    #compute correlations between fitted loadings and node centrality
    outstruct$graph_v_factor[[method]]$corr_v_fitted <- outstruct$graph_v_factor[[method]]$metric_v_loadings %>%
        group_by(graphNum, measure, factor) %>% do({
              #beware of occasional lack of variation in correlation (e.g., happens when there is simple structure in multi-factor models)
              if (sd(.$fittedloading) < 1e-8 || sd(.$value) < 1e-3) { cv = NA } else { cv = cor(.$fittedloading, .$value) }
              data.frame(graphNum=.$graphNum[1], measure=.$measure[1], cv=cv, stringsAsFactors=FALSE)
            }) %>% #summarize(cv=cor(loading,value)) %>%
        group_by(measure, factor) %>% summarize(mcv=mean(cv, na.rm=TRUE))
    
  }

  if (!saveSimsemout) { outstruct$simsemout <- NULL } #this takes up a lot of disk space, especially with updated package
  return(outstruct)
}


demomaster <- function(nexamples=2, nindicators=20, nfactors=2, nreplications=200, n=400, 
    loadings="random", loadingsampler=seq(.4, .95, .05),
    prop_residcor=0.0, mag_residcor=function() { runif(1, 0.3, 0.7) },
    prop_crossload=0.0, mag_crossload=function() { runif(1, 0.3, 0.7) },
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL, factor_correlation=0, ...) {
  
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
    } else {
      for (i in 1:nfactors) {
        lambda[i,(i-1)*perfactor + 1:perfactor] <- 0.5 #just a dummy here to make the cross-loadings code work. This will be overridden in the randomization within replications
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
            if (class(loadingsampler) == "function") {
              lambda[i,(i-1)*perfactor + 1:perfactor] <- loadingsampler(perfactor) #draw from some function/distribution
            } else {
              lambda[i,(i-1)*perfactor + 1:perfactor] <- sample(loadingsampler, perfactor, replace=TRUE)
            }
          }      
        }

        fvar <- rep(1.0, nfactors) #standardized factor variance
        psi <- matrix(factor_correlation, nrow=nfactors, ncol=nfactors) #allow for uniform factor correlation
        diag(psi) <- fvar
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
          #compute item residual variances assuming equal observed variances
          #note that we should be worried if ivar != 1.0 in a standardized solution because it can easily given negative variances
          errorvars <- computeResidvar(targetitemvar=ivar, lambda, fvar=fvar)          
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
        
        sims <- simCFAGraphs(m, nreplications=nreplications, n=n, ...)
        
        return(sims)
      })
}


##analysis functions

summarize_loadings_convergence <- function(examplestruct) {
  require(tidyverse)
  require(knitr)
  require(cowplot)
  
  ##just get a mean and SD of correlations between each nodal measure and factor loadings
  edgebasis <- c("EBICglasso", "pcor") #methods for computing graphs
  metrics <- c("strength", "closeness", "betweenness")
  retlist <- list()
  for (eb in edgebasis) {
    corrvgraph <- do.call(rbind, lapply(examplestruct, function(example) { example$graph_v_factor[[eb]]$corr_v_fitted }))
    
    #note that $metric_v_loadings has detailed statistics. In particular, it is 6000 rows per example: 200 replications x 3 metrics x 10 indicators
    corrvgraph_detailed <- do.call(rbind, lapply(1:length(examplestruct), function(ex) { 
          df <- examplestruct[[ex]]$graph_v_factor[[eb]]$metric_v_loadings %>% mutate(value=ifelse(measure=="betweenness", log10(value + 1.1) + 1, value))  #mutate(value=ifelse(measure=="betweenness", sqrt(value), value))
          #df$example <- ex
          
          #NB: if we group by factor and graph number (replication), then correlations are based on 10 samples (the 10 indicators)
          #this could be risky in terms of noisiness of the association
          #we could prefer computing correlations across all indicators and replications since these were drawn from the population model.
          dfout <- df %>% group_by(factor, measure, graphNum) %>% summarize(r=cor(fittedloading, value)) %>% summarize(m_overreps=mean(r), sd_overreps=sd(r))
          #dfout <- df %>% group_by(factor, measure) %>% summarize(r=cor(fittedloading, value)) #correlations across indicators and replications
        }))
    
    #summarize(mr=mean(r), sdr=sd(r))
    
    #summary table
    cat("\n\nAssociations between factor loadings and graph metrics for ", eb, "\n")
    #print(kable(corrvgraph_detailed %>% group_by(measure, factor) %>% summarize(m_overexamples=mean(m_overreps), sd_overexamples=sd(m_overreps),
    print(kable(corrvgraph_detailed %>% group_by(measure) %>% summarize(m_overexamples=mean(m_overreps), sd_overexamples=sd(m_overreps), #aggregate over factors
          m_wiexample_sd=mean(sd_overreps), sd_wiexamples_sd=sd(sd_overreps)), 
        digits=4, caption = paste0("Association of ", eb, " graph measures with 1-factor CFA loadings")))
    
    
    #compute plots of association between loadings and metrics
    glist <- list()
    for (m in metrics) {
      smat <- do.call(rbind, lapply(1:length(examplestruct), function(ex) {
            #brep <- filter(examplestruct[[ex]]$graph_v_factor[[eb]]$metric_v_loadings, node=="y1" & measure==m & fittedloading < 1)
            brep <- filter(examplestruct[[ex]]$graph_v_factor[[eb]]$metric_v_loadings, measure==m & fittedloading < 1)
            if (m == "betweenness") { brep$value <- log10(brep$value + 1.1) + 1 }
            brep$example <- ex
            brep
          }))
      
      bb <- smat %>% select(graphNum, value, node, example, factor, fittedloading) %>% spread(key=factor, value=fittedloading) %>%
        filter(node %in% c("y1", "y2", "y3")) #for ease of display
      cv <- bb %>% group_by(node) %>% summarize(cv=cor(value, f1)) %>% mutate(x=min(bb$f1) + 0.1*diff(range(bb$f1)), y=max(bb$value) - 0.3*diff(range(bb$value)))
      g <- ggplot(bb, aes(x=f1, y=value)) + geom_point(alpha=0.5) + stat_smooth(method="lm") + xlab("Fitted factor loading") +
        ylab(m) + theme_bw(base_size=15) + geom_text(data=cv, aes(x=x, y=y, label=paste0("r = ", round(cv, 2)))) + #annotate(geom="text", x = 0.4, y=1.1, label=paste0("r = ", round(cv, 2))) + theme_cowplot(font_size=20) +
        facet_wrap(~ node)
      glist[[m]] <- g
    }
    
    retlist[[eb]] <- list(corrvgraph_detailed = corrvgraph_detailed, corrvgraph = corrvgraph, glist=glist)
    #do.call(plot_grid, c(glist, ncol=1))
  }
  
  return(retlist)  
}



#this returns 400 x 100. 400 is the 20x20 off-factor matrix (target indicators x off-factor indicators)
#100 is the number of replications in this example
get_off_factor_dist <- function(example, method="EBICglasso", perfactor=10) {
  require(reshape2)
  
  nfactors <- dim(example$adjmats[[method]]$concat)[2] / perfactor
  
  #replications x variables x off-factor edge
  off_factor_dist <- plyr::aaply(example$adjmats[[method]]$concat, 1, function(replication) {
      nvar <- ncol(replication)
      off_factor <- matrix(NA, nrow=nvar, nvar - perfactor) #no columns for on-factor
      for (f in 1:nfactors) {
        offset <- (f-1)*perfactor
        for (o in 1:perfactor) {
          off_factor[o+offset,] <- replication[o+offset, (1:nvar)[-1*(offset + 1:perfactor)] ]
        }
      }
      return(off_factor)
    })
  
  off_sum <- apply(off_factor_dist, c(1,2), function(v) { sum(abs(v)) })
  
#  off_factor_dist_old <- apply(example$adjmats[[method]]$concat, 1, function(replication) {
#      off_factor <- matrix(NA, nrow=ncol(replication), ncol(replication) - perfactor)
#      for (o in 1:nrow(off_factor)) {
#        if (o <= 10) {
#          block <- replication[o, 11:20] #this is hard-coded for 10 indicators per factor
#        } else {
#          block <- replication[o, 1:10]
#        }
#        off_factor[o,] <- block
#      }
#      return(off_factor)
#    })
  
  #vv_off <- reshape2::melt(off_factor_dist, varnames=c("element", "replication"))
  
  vv_off <- reshape2::melt(off_sum, varnames=c("graphNum", "node"))
  
  on_factor_dist <- plyr::aaply(example$adjmats[[method]]$concat, 1, function(replication) {
      on_factor <- matrix(NA, nrow=ncol(replication), perfactor - 1)
      for (f in 1:nfactors) {
        offset <- (f-1)*perfactor
        for (o in 1:perfactor) {
          on_factor[o+offset,] <- replication[o+offset, (offset + 1:perfactor)[-1*(o)] ]
        }
      }
      return(on_factor)
    })
  
#  on_factor <- apply(example$adjmats[[method]]$concat, 1, function(replication) {
#      on_factor <- matrix(NA, nrow=ncol(replication), ncol=ncol(replication)/nfactors - 1) #always remove the self-correlation element
#      for (o in 1:length(on_factor)) {
#        if (o > 10) {
#          block <- replication[o, 11:20] #this is hard-coded for 10 indicators per factor
#        } else {
#          block <- replication[o, 1:10]
#        }
#        #wvec <- block[lower.tri(block)]
#        on_factor[o,] <- block[-o] #mean(abs(block[block != 0])) #remove zero self-correlation (diagonal) from mean
#      }
#      return(on_factor)
#    })
  
  on_sum <- apply(on_factor_dist, c(1,2), function(v) { sum(abs(v)) })

  vv_on <- reshape2::melt(on_sum, varnames=c("graphNum", "node"))
  return(list(off=vv_off, on=vv_on))
}