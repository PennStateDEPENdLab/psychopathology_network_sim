#item variance formula: 
#  var(y_i) = loading_i^2 * var(factor) + var(e_i).
#
#So, if we have loading^2 under a standardized factor model (var = 1.0), and
#  a target item variance, then we can compute the error variance as:
#  var(e_i) = var(y_i) - loadingi^2
#
#Example: factor loading of 0.9, but target item variance of 1.0
#  var(e_i) = 1.0 - .9^2 = .19
computeResidvar <- function(targetitemvar, floadings, fvar=1.0) { targetitemvar - floadings^2*fvar }

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
buildLavaanSyntax <- function(varnames, lambda, theta, psi) {
  #lambda is an nfactors x nvariables loadings matrix (lambda)
  #theta is the error/residual covariance matrix for observed dependent variables
  #psi is the covariance matrix of latent variables (factors)
  stopifnot(inherits(lambda, "matrix"))
  stopifnot(length(varnames) == nrow(lambda))
  stopifnot(nrow(psi) == ncol(lambda)) #fvar should have the same number of elements as rows of loadings
  stopifnot(nrow(psi) == ncol(psi)) #psi should be square (covariance matrix)
  syntax=c()
  fitsyntax <- c() #syntax for fitting simulated replications (doesn't get anything for free)
  #setup factor structure
  for (i in 1:ncol(lambda)) {
    syntax <- c(syntax, paste0("f", i, " =~ ", paste(lambda[,i], varnames, sep='*', collapse=" + ")))
    fitsyntax <- c(fitsyntax, paste0("f", i, " =~ NA*", paste(varnames, collapse=" + "))) #force first loading to be free by prepending NA (standardized solution)
  }
  
  #setup factor variance structure (psi)
  for (i in 1:nrow(psi)) {
    for (j in 1:ncol(psi)) {
      if (i > j) next #only lower triangle 
      if (psi[i,j] != 0) { syntax <- c(syntax, paste0("f", i, " ~~ ", psi[i,j], "*", "f", j)) } #only add syntax if non-zero
      if (i==j) { fitsyntax <- c(fitsyntax, paste0("f", i, " ~~ ", psi[i,j], "*", "f", j)) } #add standardized factor variances
    }
  }
  
  #setup error variance structure (theta)
  for (i in 1:nrow(theta)) {
    for (j in 1:ncol(theta)) {
      if (i > j) next #only lower triangle 
      if (theta[i,j] != 0) { syntax <- c(syntax, paste0(varnames[i], " ~~ ", theta[i,j], "*", varnames[j])) } #only add syntax if non-zero      
    }
  }
  return(list(simsyntax=paste(syntax, collapse="\n"), fitsyntax=paste(fitsyntax, collapse="\n")))
}

simCFAGraphs <- function(model, nreplications, n, graphmethods=c("EBICglasso", "pcor")) { #, "IsingFit"
  require(simsem)
  require(parallel)
    
  syntax <- buildLavaanSyntax(model$varnames, model$lambda, model$theta, model$psi) #build lavaan syntax for simulation
  
  #note that by default, the fitted models are provided only with a configural model on which free parameters are estimated
  #browser()
  simStruct <- simsem::sim(nRep=nreplications, model=syntax$fitsyntax, n=n, generate=syntax$simsyntax, 
      lavaanfun = "cfa", outfundata=savedat, multicore=TRUE) #TRUE)
  
  #obtain raw simulated data for network analysis
  dlist <- getExtraOutput(simStruct) #only seems to work with convergence of estimated model?
  #dlist <- simStruct@extraOut #manual (always works)
  
  outstruct <- list()
  outstruct[["specification"]] <- list(model=model, nreplications=nreplications, n=n, graphmethods=graphmethods, syntax=syntax) #snapshot of simulation structure
  outstruct[["simsemout"]] <- simStruct
  outstruct[["simdata"]] <- dlist
  
  for (method in graphmethods) {
    glist <- mclapply(dlist, function(df) {
          estimateNetwork(df, default=method)
        })
    
    outstruct[["graphs"]][[method]] <- glist #save raw graphs into output
    
    #compute centrality metrics (just the handful supported in bootnet
    Long <- mclapply(1:length(glist), function(g) { df <- centralityTable(glist[[g]], standardized=TRUE); df$graphNum <- g; return(df) })
    LongAll <- do.call(rbind, Long) #includes metrics for all replications (for graphing means and variation/uncertainty
    
    outstruct[["gmetrics"]][[method]] <- LongAll #save centrality metrics into output
  }
  
  #tests of different models
  #partial correlation
  #partcorr <- dlist[[1]] %>% qgraph::cor_auto() %>% corpcor::cor2pcor()
  
  #binarized ising model
  #ising <- dlist[[1]] %>% bootnet::binarize() %>% IsingFit::IsingFit()
  
  #dlist[[1]] %>% qgraph::cor_auto() %>% qgraph::EBICglasso()
    
  return(outstruct)
}

##First demonstration: correlation of one-factor CFA loadings with centrality measures
##a function to simulate multiple one-factor CFAs:
## 10, 15, or 20 indicators
## a basis for sampling loadings like uniform from seq(.4, .95, .05)
## a number of replications for each setting.

demo1 <- function(nexamples=10, nindicators=10, nreplications=200, n=400, loadingsampler=seq(.4, .95, .05), errorvars="eqresidvar", ivar=0.5) {  
  lapply(1:nexamples, function(i) {
        lambda <- matrix(rbind(sample(loadingsampler, nindicators, replace=TRUE)))
        varnames <- paste0("y", 1:nrow(lambda))
        fvar <- 1.0 #factor variance (standardized)
        if (errorvars=="eqresidvar") {
          errorvars <- rep(ivar, length(varnames)) #fixed error specification          
        } else if (errorvars=="eqvars") {
          errorvars <- computeResidvar(targetitemvar=ivar, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances          
        }
        
        theta <- diag(as.vector(errorvars)) #0 resid cov by default
        rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work
        
        #populate a few error correlations
        #theta <- addErrorCor(theta, c("y1", "y2"), 0.5)
        #theta <- addErrorCor(theta, c("y1", "y4"), 0.6)
        
        psi <- diag(fvar) #zero covariance in factor structure at the moment
        
        #model specification structure. Currently just wrapping up individual arguments above into a list
        m <- list(
            varnames=varnames, #vector of variable names
            lambda=lambda, #nitems x nfactors loadings matrix
            theta=theta, #covariance (residual) matrix for observed items
            psi=psi #covariance matrix for factors
        )
        
        sims <- simCFAGraphs(m, nreplications=nreplications, n=n)
        
        loadmaster <- data.frame(poploading=as.vector(lambda), node=varnames)
        
        #I should come back and populate other methods
        sims$gmetrics$EBICglasso <- join(sims$gmetrics$EBICglasso, loadmaster, by="node")
        
        fitted_loadings <- select(sims[["simsemout"]]@coef, starts_with("f1=~"))
        fitted_loadings$graphNum <- 1:nrow(fitted_loadings)
        names(fitted_loadings) <- sub("f1=~", "", names(fitted_loadings), fixed=TRUE)
        mm <- gather(fitted_loadings, key="node", value="fittedloading", -graphNum)
        
        sims$gmetrics$EBICglasso <- join(sims$gmetrics$EBICglasso, mm, by=c("node", "graphNum"))
        
        return(sims)
      })
}



#LEFTOVERS
#cfa1 <- '
#    f1 =~ 1*y1 + 0.6*y2 + 0.7*y3 + 0.8*y4 + 0.9*y5 + 0.6*y6 + 0.45*y7 + 0.9*y8 + 1*y9
#    
#    f1 ~~ 1.0*f1 #1.0 factor variance (standardize)
#    
#    #similar residual variances for indicators
#    #item variances are loading x  
#    
#    y1 ~~ 0.7*y1
#    y2 ~~ 0.9*y2
#    y3 ~~ 0.8*y3
#    y4 ~~ 0.7*y4
#    y5 ~~ 0.6*y5
#    y6 ~~ 0.95*y6
#    y7 ~~ 0.6*y7
#    y8 ~~ 0.7*y8
#    y9 ~~ 0.8*y9
#    '
#
#analyzeModel <- '
#    f1 =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9
#    '