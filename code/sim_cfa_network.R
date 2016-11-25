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

##demonstration 1: centrality is inversely proportional to factor loadings 
#example one-factor CFA
lambda <- as.matrix(rbind(c(0.7, 0.6, 0.7, 0.8, 0.9, 0.6, 0.5, 0.9, 0.4)))
varnames <- paste0("y", 1:ncol(lambda))
fvar <- 1.0 #factor variance (standardized)
errorvars <- rep(0.5, length(varnames)) #fixed error specification
#errorvars <- computeResidvar(targetitemvar=1.0, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances 
theta <- diag(as.vector(errorvars)) #0 resid cov initially
rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work

#populate a few error correlations
#theta <- addErrorCor(theta, c("y1", "y2"), 0.5)
#theta <- addErrorCor(theta, c("y1", "y4"), 0.6)

psi <- diag(fvar) #zero covariance in factor structure at the moment
dimnames(psi) <- list(f=paste0("f", 1:nrow(psi)), f=paste0("f", 1:ncol(psi)))

#syntax <- buildLavaanSyntax(varnames, loadings, theta, psi)

#model specification structure. Currently just wrapping up individual arguments above into a list
model2 <- list(
    varnames=varnames, #vector of variable names
    lambda=lambda, #nitems x nfactors loadings matrix
    theta=theta, #covariance (residual) matrix for observed items
    psi=psi #covariance matrix for factors
)

sims1 <- simCFAGraphs(model2, nreplications=20, n=200)

sims1$graph_v_factor$EBICglasso$corr_v_fitted
str(sims1$graph_v_factor$EBICglasso$metric_v_loadings)

#look at modularity with a one-factor CFA
g <- igraph::graph_from_adjacency_matrix(sims1$adjmats[["EBICglasso"]]$concat[1,,], mode="undirected", weighted=TRUE, diag=FALSE)
g <- delete.edges(g, which(E(g)$weight < 0.01))
cluster_louvain(g)
cluster_optimal(g)
cluster_leading_eigen(g)
cluster_spinglass(g)
cluster_walktrap(g)
cluster_label_prop(g)
cluster_edge_betweenness(g)
cluster_fast_greedy(g)
largest_cliques(g)

decompose.graph(g,mode="strong")
g2 <- add.edges(g, c(1,10), weight=.1)





#start with EBICGLASSO
#need to extract centrality measures across replications and correlate them with 
str(sims$gmetrics$EBICglasso)

loadmaster <- data.frame(poploading=as.vector(lambda), node=varnames)
df <- inner_join(sims$gmetrics$EBICglasso, loadmaster, by="node")

library(lme4)
test <- lmer(value ~ poploading + (1 | graphNum), filter(df, measure=="Betweenness"))
summary(test2 <- lmer(value ~ poploading + (1 | graphNum), filter(df, measure=="Closeness")))
summary(test3 <- lmer(value ~ poploading + (1 | graphNum), filter(df, measure=="Strength")))

gmetricsAgg <- sims$gmetrics$EBICglasso %>% select(-graph, -type) %>% group_by(measure, node) %>% 
    summarize(mval=mean(value)) #%>% join(loadmaster, by="node")

df2 <- merge(gmetricsAgg, loadmaster, by="node")
df2 %>% group_by(measure) %>% summarize(cor(mval, poploading))

##so when we aggregate correlate the 9 population loadings, we get huge correlations
##when we maintain equal item variances, the corrs for closeness and strength are around ~.95
##when we set equal residual variances, the corrs are essentially 1
#> df2 %>% group_by(measure) %>% summarize(cor(mval, poploading))
# A tibble: 3 × 2
#       measure `cor(mval, poploading)`
#       <fctr>                   <dbl>
#1 Betweenness               0.9509512
#2   Closeness               0.9993568
#3    Strength               0.9996537

#what about correlation between fitted factor loadings and centrality
#fitted loadings: sims[["simsemout"]]@coef)
floadings <- select(sims[["simsemout"]]@coef, starts_with("f1=~"))
floadings$graphNum <- 1:nrow(floadings)
names(floadings) <- sub("f1=~", "", names(floadings), fixed=TRUE)
mm <- gather(floadings, key="node", value="poploading", -graphNum)

df3 <- inner_join(sims$gmetrics$EBICglasso, mm, by=c("node", "graphNum"))

df3 %>% group_by(measure) %>% summarize(cor(value, poploading))




#do things blow up when we allow for big correlations?
thetacorlist <- list(
    list("y1", "y2", 0.8),
    list("y1", "y3", 0.1),
    list("y4", "y6", 0.9),
    list("y7", "y9", 0.7)
)
dd3 <- demo1(thetacorlist=thetacorlist)

qgraph(dd3[[1]]$adjmats$EBICglasso$average)
qgraph(dd3[[1]]$adjmats$pcor$average)
qgraph(dd3[[1]]$adjmats$pearson$average)

#yes, this yields a steep departure from the fitted loadings, even if we allow the residual variances into the estimation
dd3[[1]]$gmetrics$EBICglasso %>% group_by(graphNum, measure) %>% summarize(cv=cor(fittedloading,value)) %>% 
    group_by(measure) %>% summarize(mcv=mean(cv))

#but what if we revert to a pearson correlation basis? shouldn't be much better in principle since there are
#effectively two causes

dd3[[1]]$gmetrics$EBICglasso %>% group_by(graphNum, measure) %>% summarize(cv=cor(fittedloading,value)) %>% 
    group_by(measure) %>% summarize(mcv=mean(cv))


#what if we allowed for the residual loadings in the CFA. Would the factor loadings come back into range?
#doesn't seem to improve much


