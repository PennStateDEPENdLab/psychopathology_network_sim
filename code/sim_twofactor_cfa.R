#purpose: prototype large two-factor CFA simulation that is analyzed using network methods.
setwd(file.path(getMainDir(), "psychopathology_network_sim"))
source("code/psychopathology_sim_functions.R")
library(simsem)
library(bootnet)
library(qgraph)
library(dplyr)
library(ggplot2)
library(abind)
library(doParallel)

#####
#two-factor setup
lambda <- as.matrix(rbind(
        c(0.7, 0.6, 0.7, 0.8, 0.9, 0.6, 0.45, 0.9, 0.7, rep(0, 9)),
        c(rep(0, 8), 0.6, 0.7, 0.6, 0.7, 0.8, 0.9, 0.6, 0.45, 0.9, 0.7)))

varnames <- paste0("y", 1:ncol(lambda))

fvar <- c(1.0, 1.0) #factor variance (standardized)
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

sims_2fac <- simCFAGraphs(model1, nreplications=100, n=200)

summary(sims_2fac$simsemout)
#df <- sims$simdata[[1]]
#save(file="lvsimdata_forclass.RData", df)

pdf("figures/est 2fac ebic avg network error cor.pdf", width=8, height=8)
#just plot the average adjacency directly
#qgraph(sims_2fac$adjmats$EBICglasso$average)
qgraph(sims_2fac$adjmats$EBICglasso$concat[2,,])
dev.off()

summary(sims_2fac$simsemout)
sims_2fac$graph_v_factor$EBICglasso$corr_v_fitted
str(sims_2fac$graph_v_factor$EBICglasso$metric_v_loadings)
table(sims_2fac$graph_v_factor$EBICglasso$metric_v_loadings$poploading)

sims_2fac$graph_v_factor$EBICglasso$corr_v_fitted
#basically need to switch over to igraph to get a feel for structure (subgraphs)

##closeness is at 0 in some cases
g <- igraph::graph_from_adjacency_matrix(sims_2fac$adjmats[["EBICglasso"]]$concat[33,,], mode="undirected", weighted=TRUE, diag=FALSE)
g <- delete.edges(g, which(E(g)$weight < 0.01))
cluster_louvain(g)
cluster_optimal(g)
g2 <- add.edges(g, c(1,10), weight=.1)
cluster_leading_eigen(g2)

#this is probably the most sane: a connected component has edges between all vertices (fully connected)
#the strong versus weak does not apply here -- that's just for directed graphs
dec <- decompose.graph(g, min.vertices=3)
strength1 <- strength(dec[[1]])
cluster_louvain(dec[[1]])

#contrived example
g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#g <- add_edges(g, c(1,6, 1,11, 6, 11))
wtc <- cluster_walktrap(g)

sims_2fac$graphs$EBICglasso[[33]]$graph
centrality_auto(sims_2fac$graphs$EBICglasso[[33]])


igraph::closeness(g, normalized=FALSE, weights=1/igraph::E(g)$weight)
igraph::betweenness(g, normalized=FALSE, weights=1/igraph::E(g)$weight)

## HYPOTHESES:
## 1) In 2-factor structure, betweenness will moderately capture cross-loading
## 2) In 2-factor structure, the magnitude of a residual correlation will strongly correlate with the betweenness of the two nodes connected.
## 3) In n-factor structure, strength or degree should remain equally associated (as a 1-fac) with factor loading because
##     the edges within a factor reflect the relevant factor loadings and there will not be edges to other factors (absent cross-loadings).
## 4) In 2-factor structure, if we decompose the graph into strongly connected subcomponents (either through community algorithm or through
##     a largest connected component decomposition), then all of the 1-factor results will hold (metrics within factors).

## TESTS:
## 1) Scale up cross loading in a 2-factor structure that otherwise remains static
## 2) Look at strength and degree in 2-, 3-, and 4-factor models. For degree, perhaps threshold at 0.4 density -- or even use a range?
## 3) For the decomposition analysis, use modularity analysis or largest connected component decomposition, then compute metrics within subgraphs

#first up: cross-loadings
x <- demomaster(nexamples=2, nindicators=20, nfactors=2, nreplications=200, n=400, 
    loadings=0.9, prop_crossload=0.2, mag_crossload=function() { runif(1, 0.3, 0.7) },
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)
  
firstex <- x[[1]]
str(firstex, max=1)
summary(firstex$simsemout)

firstex$graph_v_factor$EBICglasso$corr_v_fitted
cat(firstex$specification$syntax$fitsyntax)
cat(firstex$specification$syntax$simsyntax)


x3 <- demomaster(nexamples=1, nindicators=20, nfactors=2, nreplications=200, n=400, 
    loadingsampler=seq(.4, .95, .05), errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

#basic setup: correlations of graph with factor for 1-4 factors
x1 <- demomaster(nexamples=1, nindicators=20, nfactors=1, nreplications=200, n=400, 
    loadingsampler=seq(.4, .95, .05), errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

x2 <- demomaster(nexamples=1, nindicators=20, nfactors=2, nreplications=200, n=400, 
    loadingsampler=seq(.4, .95, .05), errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

x3 <- demomaster(nexamples=1, nindicators=20, nfactors=3, nreplications=200, n=400, 
    loadingsampler=seq(.4, .95, .05), errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

x4 <- demomaster(nexamples=1, nindicators=20, nfactors=4, nreplications=200, n=400, 
    loadingsampler=seq(.4, .95, .05), errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

#hypothesis 3 corroborated
x1[[1]]$graph_v_factor$EBICglasso$corr_v_fitted
#summary(x[[1]]$simsemout)
x2[[1]]$graph_v_factor$EBICglasso$corr_v_fitted
x3[[1]]$graph_v_factor$EBICglasso$corr_v_fitted
x4[[1]]$graph_v_factor$EBICglasso$corr_v_fitted

##
## 1) In 2-factor structure, betweenness will capture cross-loading

##need to lock in a given factor structure, then vary the cross-loading magnitude.

#####
#consistent 2-factor structure
lambda <- as.matrix(rbind(
        c(rep(0.6, 9), rep(0, 9)), #f1
        c(rep(0, 9), rep(0.6, 9)))) #f2

varnames <- paste0("y", 1:ncol(lambda))

fvar <- c(1.0, 1.0) #factor variance (standardized)
errorvars <- rep(0.5, ncol(lambda)) #fixed error specification
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

#add one cross-loading...
mlist <- lapply(1:500, function(i) {
      mm <- model1
      mm$lambda[2, 2] <- runif(1, 0.2, 0.8) #add cross-loading of y2 onto f2
      sims <- simCFAGraphs(mm, nreplications=1, n=200, thetastart=TRUE, parallel=0)      
      bmat <- filter(sims$graph_v_factor$EBICglasso$metric_v_loadings, node=="y2" & measure=="betweenness")
      #mload <- cor.test(~ fittedloading + value, filter(bmat, factor=="f1"))$estimate #correlation of node betweenness and core factor
      #cload <- cor.test(~ fittedloading + value, filter(bmat, factor=="f2"))$estimate #correlation of node betweenness and core factor
      bb <- bmat %>% select(graphNum, value, factor, fittedloading) %>% spread(key=factor, value=fittedloading) #graphNum contains the identifying column
      #summary(lm(value ~ f1 + f2, bb))
      #return(c(mload=mload, cload=cload))
      return(select(bb, -graphNum))
    })

mdf <- do.call(rbind, mlist)
cor(mdf$value, mdf$f1)
cor(mdf$value, mdf$f2) #weak, but significant, support

summary(lm(value ~ f1 + f2, mdf))

sims <- simCFAGraphs(model1, nreplications=100, n=200, thetastart=TRUE)

summary(sims$simsemout)

#more direct is residual correlation

cl <- makeCluster(8) #defaults to PSOCK cluster
#clusterExport(cl, c("estimateNetwork"))#, "cd4.mle"))
clusterExport(cl, c("simCFAGraphs", "model1", "addErrorCor", "buildLavaanSyntax", "savedat"))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(qgraph))
clusterEvalQ(cl, library(igraph))
clusterEvalQ(cl, library(tidyr))


mlist <- parLapply(cl, 1:500, function(i) {
      mm <- model1
      mm$theta <- addErrorCor(mm$theta, c("y1", "y10"), runif(1, 0.2, 0.8))
      
      sims <- simCFAGraphs(mm, nreplications=1, n=200, thetastart=TRUE, parallel=0)
      
      bmat <- filter(sims$gmetrics$EBICglasso, node %in% c("y1", "y10") & measure=="betweenness") #need to grab magnitude of residual correlation...
      bb <- bmat %>% select(graphNum, node, value, rvar, rcorr_fitted) %>% spread(key=rvar, value=rcorr_fitted) #graphNum contains the identifying column
      #summary(lm(value ~ f1 + f2, bb))
      #return(c(mload=mload, cload=cload))
      return(select(bb, -graphNum))
    })

stopCluster(cl)

#should probably reshape to allow for MLM where nodes are random
#initial corroboration for hypothesis 2
mm <- do.call(rbind, mlist)
my1 <- filter(mm, node=="y1")
summary(lm(value ~ y10, my1)) #betweenness of y1 predicted by residual correlation with y10
cor.test(~ value + y10, my1) #r ~ 0.55

my10 <- filter(mm, node=="y10")
summary(lm(value ~ y1, my10)) #betweenness of y10 predicted by residual correlation with y1
cor.test(~ value + y1, my10)

#x2 <- demomaster(nexamples=1, nindicators=20, nfactors=2, nreplications=200, n=400, 
#    loadings="random", loadingsampler=seq(.4, .95, .05), prop_crossload=0.2, mag_crossload=function() { runif(1, 0.3, 0.7) },
#    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)
#
#x2[[1]]$graph_v_factor$EBICglasso$corr_v_fitted
#summary(x2[[1]]$simsemout)
