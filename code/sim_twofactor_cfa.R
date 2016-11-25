#purpose: prototype large two-factor CFA simulation that is analyzed using network methods.
setwd(file.path(getMainDir(), "psychopathology_network_sim"))
source("code/psychopathology_sim_functions.R")
library(simsem)
library(bootnet)
library(qgraph)
library(dplyr)
library(ggplot2)
library(abind)

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
## 1) In 2-factor structure, betweenness will capture cross-loading
## 2) In 2-factor structure, the magnitude of a residual correlation will correlate with the betweenness of the two nodes connected.
## 2) In n-factor structure, strength or degree should remain equally associated (as a 1-fac) with factor loading because
##     the edges within a factor reflect the relevant factor loadings and there will not be edges to other factors (absent cross-loadings).
## 3) In 2-factor structure, if we decompose the graph into strongly connected subcomponents (either through community algorithm or through
##     a largest connected component decomposition), then all of the 1-factor results will hold (metrics within factors).

## TESTS:
## 1) Scale up cross loading in a 2-factor structure that otherwise remains static
## 2) Look at strength and degree in 2-, 3-, and 4-factor models. For degree, perhaps threshold at 0.4 density -- or even use a range?
## 3) For the decomposition analysis, use modularity analysis or largest connected component decomposition, then compute metrics within subgraphs

#first up: cross-loadings
x <- democrossload(nexamples=2, nindicators=20, nfactors=2, nreplications=200, n=400, 
    loadings=0.9, prop_crossload=0.2, mag_crossload=function() { runif(1, 0.3, 0.7) },
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)
  
firstex <- x[[1]]
str(firstex, max=1)
summary(firstex$simsemout)

firstex$graph_v_factor$EBICglasso$corr_v_fitted
cat(firstex$specification$syntax$fitsyntax)
cat(firstex$specification$syntax$simsyntax)


x <- demomaster(nexamples=1, nindicators=20, nfactors=4, nreplications=20, n=400, 
    loadings=0.9, prop_crossload=0.2, mag_crossload=function() { runif(1, 0.3, 0.5) },
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

summary(x[[1]]$simsemout)