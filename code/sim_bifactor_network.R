#purpose: prototype large one-factor CFA simulation that is analyzed using network methods.
#on the network side, largely follow the Epskamp bootnet paper
setwd(file.path(getMainDir(), "psychopathology_network_sim"))
source("code/psychopathology_sim_functions.R")
library(simsem)
library(bootnet)
library(qgraph)
library(dplyr)
library(ggplot2)
library(abind)
library(tidyr)

#####
#bifactor setup
lambda <- as.matrix(rbind(
        rep(0.6, 18), #general factor
        c(rep(0.6, 9), rep(0, 9)), #specific 1
        c(rep(0, 9), rep(0.6, 9)))) #specific 2
        #c(0.7, 0.6, 0.7, 0.8, 0.8, 0.6, 0.45, 0.8, 0.7, rep(0, 9)), #specific 1
        #c(rep(0, 9), 0.7, 0.6, 0.7, 0.8, 0.8, 0.6, 0.45, 0.8, 0.7))) #specific 2

varnames <- paste0("y", 1:ncol(lambda))

fvar <- c(1.0, 1.0, 1.0) #factor variance (standardized)
#errorvars <- rep(0.5, length(loadings)) #fixed error specification
errorvars <- computeResidvar(targetitemvar=1.0, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances 
theta <- diag(as.vector(errorvars)) #0 resid cov initially
rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work

#populate a few error correlations
#theta <- addErrorCor(theta, c("y1", "y2"), 0.5)
#theta <- addErrorCor(theta, c("y1", "y4"), 0.6)

psi <- diag(fvar) #zero covariance in factor structure at the moment
dimnames(psi) <- list(f=paste0("f", 1:nrow(psi)), f=paste0("f", 1:ncol(psi)))
#syntax <- buildLavaanSyntax(varnames, loadings, theta, psi)

#model specification structure. Currently just wrapping up individual arguments above into a list
model1 <- list(
    varnames=varnames, #vector of variable names
    lambda=lambda, #nitems x nfactors loadings matrix
    theta=theta, #covariance (residual) matrix for observed items
    psi=psi #covariance matrix for factors
)

sims <- simCFAGraphs(model1, nreplications=100, n=200, thetastart=TRUE)

summary(sims$simsemout)

pdf("figures/est bifactor ebic avg network.pdf", width=8, height=8)
#just plot the average adjacency directly
qgraph(sims$adjmats$EBICglasso$average)
#qgraph(sims$adjmats$EBICglasso$concat[1,,])
dev.off()

#so we don't really see correlations in a bifactor structure
sims$graph_v_factor$EBICglasso$corr_v_fitted
head(filter(sims$graph_v_factor$EBICglasso$metric_v_loadings, loadingtype=="fittedloading"), n=100)

#but could the 

