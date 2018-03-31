#sim 3: Here, instead of varying the magnitude of one cross-loading while the other remains fixed, we simulate data along a 2-D grid in which
# the factor loading of y2 on f1 varies from 0.2--0.8 and its factor loading on f2 also varies from 0.2--0.8.
# The hypothesis is that centrality measures will be a joint function of the loading magnitudes such that high loadings 
# on both will be associated with high centrality. This is a more powerful demonstration of the point above that network metrics do not allow one to unmix multiple causes.


#unlike sim2, we don't force the variance explained to be equal anywhere
#though in general, the basic structure is similar, but here, we make the error variance
#0.3 throughout

library(qgraph)
library(igraph)
library(lavaan)
library(viridis)
library(dplyr)
library(lme4)

setwd("~/Data_Analysis/psychopathology_network_sim")
source("code/psychopathology_sim_functions.R")

#consistent 2-factor structure
lambda <- as.matrix(rbind(
    c(rep(0.8, 9), rep(0, 9)), #f1
    c(rep(0, 9), rep(0.8, 9)))) #f2

varnames <- paste0("y", 1:ncol(lambda))

fvar <- c(1.0, 1.0) #factor variance (standardized)
errorvars <- rep(0.3, ncol(lambda)) #fixed error specification
#errorvars <- computeResidvar(targetitemvar=1.0, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances 
theta <- diag(as.vector(errorvars)) #0 resid cov initially
rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work

psi <- diag(fvar) #zero covariance in factor structure at the moment
dimnames(psi) <- list(f=paste0("f", 1:nrow(psi)), f=paste0("f", 1:ncol(psi)))

#model specification structure. Currently just wrapping up individual arguments above into a list
model1 <- list(
  varnames=varnames, #vector of variable names
  lambda=lambda, #nitems x nfactors loadings matrix
  theta=theta, #covariance (residual) matrix for observed items
  psi=psi #covariance matrix for factors
)

cl <- makeCluster(12) #defaults to PSOCK cluster
clusterExport(cl, c("simCFAGraphs"))

registerDoParallel(cl)
clusterEvalQ(cl, library(dplyr)) #load dplyr in all workers
loading_grid <- expand.grid(f1=seq(0.2, 0.8, .05), f2=seq(0.2, 0.8, 0.05))
mlist_crossload <- foreach(i=iter(loading_grid, by='row')) %dopar% {
  mm <- model1
  mm$lambda[1, 2] <- i$f1 #runif(1, 0.2, 0.8) #random loading of y2 onto f1
  mm$lambda[2, 2] <- i$f2  #runif(1, 0.2, 0.8) #random loading of y2 onto f2
  sims <- simCFAGraphs(mm, nreplications=100, n=400, thetastart=TRUE, parallel=0, , saveLavObj = FALSE)
  return(sims)
}
stopCluster(cl)

save(mlist_crossload, file=file.path("data", "sim3_crossload_mlist.RData"))

