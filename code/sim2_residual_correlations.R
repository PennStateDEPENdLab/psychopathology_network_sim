#simulation 2 for hallquist, wright, molenaar
#take-homes from sim 1: 
#  1) strength is largely redundant with factor loadings
#  2) closeness and betweenness are extremely sensitive to spurious sampling variability with other factors
#
# simulation 2 picks up from this point and asks, what happens if we specifically manipulate residual correlations apart from factor loadings
# this relaxes the conditional independence assumption of the standard factor model, which is often a source of criticism among network modelers.

library(qgraph)
library(igraph)
library(lavaan)
library(viridis)
library(dplyr)
library(lme4)

setwd("~/Data_Analysis/psychopathology_network_sim")
source("code/psychopathology_sim_functions.R")

cl <- makeCluster(8) #defaults to PSOCK cluster
clusterExport(cl, c("simCFAGraphs"))

registerDoParallel(cl)
clusterEvalQ(cl, library(dplyr)) #load dplyr in all workers

#Simulate residual correlations on 2-D grid of cross-loadings for four items: y1, y2 (f1); y10, y11 (f2)
#Residual correlations range in strength from 0.0 - 0.8 in .05 increments while holding constant the primary factor loadings at 0.6.
#loading_grid <- expand.grid(cl1=seq(0.0, 0.8, .05), cl2=seq(0.0, 0.8, 0.05))

#on further reflection, what we really want is to pit factor loading against residual correlation
#in the later simulation, this is achieved directly by scaling holding constant variance explained for an item,
#but divvying between factor and 'u' factor (two-item latent). This gives a clean separation in terms of variance partitioning.

#what if we use an R2 of 64%?
r2 <- .64
perfac <- 10

lambda <- as.matrix(rbind(
    c(rep(sqrt(r2), perfac), rep(0, perfac)), #f1
    c(rep(0, perfac), rep(sqrt(r2), perfac)),  #f2
    c(rep(0, perfac*2)))) #u (aka f3)

lambdaconstraint <- matrix(0, nrow=nrow(lambda), ncol=ncol(lambda))
lambdaconstraint[3, c(2,11)] <- 1 #constrain these loadings to equality in fitting (loading matrix contains integer identifiers for equated params)
varnames <- paste0("y", 1:ncol(lambda))

fvar <- c(1.0, 1.0, 1.0) #factor variance (standardized)
#errorvars <- rep(0.3, ncol(lambda)) #fixed error specification
errorvars <- computeResidvar(targetitemvar=1.0, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances 
theta <- diag(as.vector(errorvars)) #0 resid cov initially
rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work

psi <- diag(fvar) #zero covariance in factor structure at the moment
dimnames(psi) <- list(f=paste0("f", 1:nrow(psi)), f=paste0("f", 1:ncol(psi)))

#model specification structure. Currently just wrapping up individual arguments above into a list
model1 <- list(
  varnames=varnames, #vector of variable names
  lambda=lambda, #nitems x nfactors loadings matrix
  lambdaconstraint=lambdaconstraint, #constraint matrix for factor loadings
  theta=theta, #covariance (residual) matrix for observed items
  psi=psi #covariance matrix for factors
)

#f1loading <- seq(0.15, 0.80, .01)
#f2loading <- sqrt(r2 - f1loading^2 + .Machine$double.eps) #to avoid slightly negative numbers
#f1loading^2 + f2loading^2

#rework the loading grid to be in step size of variance explained .01 increment
f_loading <- sqrt(seq(0, r2, .01))
u_loading <- sqrt(r2 - f_loading^2 + .Machine$double.eps)
f_loading^2 + u_loading^2
loading_grid <- data.frame(f_loading, u_loading)

#update: for identification, equate the loading of y2 and y11 on u (2-variable factor)
mlist_uniq <- foreach(i=iter(loading_grid, by='row')) %dopar% {
  mm <- model1
  
  mm$lambda[1, 2] <- i$f_loading #loading of y2 on f1
  mm$lambda[3, 2] <- i$u_loading #loading of y2 on u
  mm$lambda[2, 11] <- i$f_loading #loading of y11 on f2
  mm$lambda[3, 11] <- i$u_loading #loading of y11 on u
  sims <- simCFAGraphs(mm, nreplications=10, n=400, thetastart=TRUE, parallel=0, saveLavObj = FALSE)
  
  return(sims)
}

stopCluster(cl)

save(mlist_uniq, file=file.path("data", "u_vs_g_mlist.RData"))

#older approach where loadings are fixed, but residual correlation is random
#difficulty is that loadings have no real variation, so we can't see how they tradeoff with covariances
#perfactor <- 9
#
#lambda <- as.matrix(rbind(
#    c(rep(0.7, perfactor), rep(0, perfactor)), #f1
#    c(rep(0, perfactor), rep(0.7, perfactor)))) #f2
#
#varnames <- paste0("y", 1:ncol(lambda))
#
#fvar <- c(1.0, 1.0) #factor variance (standardized)
#errorvars <- rep(0.3, ncol(lambda)) #fixed error specification
##errorvars <- computeResidvar(targetitemvar=1.0, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances 
#theta <- diag(as.vector(errorvars)) #0 resid cov initially
#rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work
#
#psi <- diag(fvar) #zero covariance in factor structure at the moment
#dimnames(psi) <- list(f=paste0("f", 1:nrow(psi)), f=paste0("f", 1:ncol(psi)))
#
##model specification structure. Currently just wrapping up individual arguments above into a list
#model1 <- list(
#  varnames=varnames, #vector of variable names
#  lambda=lambda, #nitems x nfactors loadings matrix
#  theta=theta, #covariance (residual) matrix for observed items
#  psi=psi #covariance matrix for factors
#)


#mlist_crossload <- foreach(i=iter(loading_grid, by='row')) %dopar% {
#  mm <- model1
#  mm$theta <- addErrorCor(mm$theta, c("y1", "y10"), i$cl1)
#  mm$theta <- addErrorCor(mm$theta, c("y2", "y11"), i$cl2)
#  sims <- simCFAGraphs(mm, nreplications=5, n=400, thetastart=TRUE, parallel=0)
#  return(sims)
#}
#stopCluster(cl)