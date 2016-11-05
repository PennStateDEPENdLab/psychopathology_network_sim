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

#example one-factor CFA
lambda <- matrix(rbind(c(0.7, 0.6, 0.7, 0.8, 0.9, 0.6, 0.45, 0.9, 0.7))) #can be extended to be a 2-d matrix (factors x indicators
varnames <- paste0("y", 1:nrow(lambda))
fvar <- 1.0 #factor variance (standardized)
#errorvars <- rep(0.5, length(loadings)) #fixed error specification
errorvars <- computeResidvar(targetitemvar=1.0, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances 
theta <- diag(as.vector(errorvars)) #0 resid cov initially
rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work

#populate a few error correlations
theta <- addErrorCor(theta, c("y1", "y2"), 0.5)
theta <- addErrorCor(theta, c("y1", "y4"), 0.6)

psi <- diag(fvar) #zero covariance in factor structure at the moment

#syntax <- buildLavaanSyntax(varnames, loadings, theta, psi)

#model specification structure. Currently just wrapping up individual arguments above into a list
model1 <- list(
    varnames=varnames, #vector of variable names
    lambda=lambda, #nitems x nfactors loadings matrix
    theta=theta, #covariance (residual) matrix for observed items
    psi=psi #covariance matrix for factors
)
       
sims <- simCFAGraphs(model1, nreplications=100, n=200)

  ff <- sims$specification$syntax$fitsyntax
  ff <- paste(ff, "y1 ~~ y2", sep="\n")
  ff <- paste(ff, "y1 ~~ y4", sep="\n")
  
  f1 <- cfa(ff, data=sims$simdata[[1]])
  summary(f1, standardized=TRUE)

#trying to get a picture of the model
#popmodel <- cfa(model=analyzeModel, data=dlist[[1]])
#popmodel <- cfa(model=cfa1)
#pdf("pop model.pdf", width=10, heigh=8)
#semPaths(popmodel,  "std", curvePivot = TRUE, edge.label.cex = 0.65, fixedStyle = c("black",1))
#dev.off()

#semPaths(simStruct@pt)
#summary(simStruct)

#one pcor graph (first replication)
pcorNet <- sims$graphs$pcor[[1]] 
pdf("figures/est network.pdf", width=8, height=8)
plot(pcorNet, layout="spring", labels=TRUE)
dev.off()

#EBICglasso centrality plots
pdf("figures/GLASSO Centrality Density.pdf", width=6, height=15)
ggplot(sims$gmetrics$EBICglasso, aes(x=value)) + geom_density() + facet_grid(node ~ measure, scales="free") + theme_bw(base_size=20)
dev.off()

#summary stats
gsums <- sims$gmetrics$EBICglasso %>% group_by(measure, node) %>% summarize(m=mean(value), sd=sd(value), se=plotrix::std.error(value))

pdf("figures/GLASSO Centrality Summaries.pdf", width=15, height=6)
ggplot(gsums, aes(x=node, y=m, ymin=m-se, ymax=m+se)) + facet_wrap(~measure, scales="free_y") + geom_pointrange() + 
    ylab("Mean of graph metric") + xlab("Indicator/Node") + theme_bw(base_size=20)
dev.off()

#average adjacency matrix (simple correlation)
allAdj <- do.call(abind, list(lapply(sims$simdata, cor), along=0))
avgAdj <- apply(allAdj, c(2,3), mean)

##demonstration 1: centrality is inversely proportional to factor loadings 
#example one-factor CFA
lambda <- matrix(rbind(c(0.7, 0.6, 0.7, 0.8, 0.9, 0.6, 0.5, 0.9, 0.4)))
varnames <- paste0("y", 1:nrow(lambda))
fvar <- 1.0 #factor variance (standardized)
errorvars <- rep(0.5, length(varnames)) #fixed error specification
#errorvars <- computeResidvar(targetitemvar=1.0, lambda, fvar=fvar) #compute item residual variances assuming equal observed variances 
theta <- diag(as.vector(errorvars)) #0 resid cov initially
rownames(theta) <- colnames(theta) <- varnames #necessary for addErrorCor to work

#populate a few error correlations
#theta <- addErrorCor(theta, c("y1", "y2"), 0.5)
#theta <- addErrorCor(theta, c("y1", "y4"), 0.6)

psi <- diag(fvar) #zero covariance in factor structure at the moment

#syntax <- buildLavaanSyntax(varnames, loadings, theta, psi)

#model specification structure. Currently just wrapping up individual arguments above into a list
model2 <- list(
    varnames=varnames, #vector of variable names
    lambda=lambda, #nitems x nfactors loadings matrix
    theta=theta, #covariance (residual) matrix for observed items
    psi=psi #covariance matrix for factors
)

sims <- simCFAGraphs(model2, nreplications=200, n=200)

#start with EBICGLASSO
#need to extract centrality measures across replications and correlate them with 
str(sims$gmetrics$EBICglasso)

loadmaster <- data.frame(poploading=as.vector(lambda), node=varnames)
df <- join(sims$gmetrics$EBICglasso, loadmaster, by="node")

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

df3 <- join(sims$gmetrics$EBICglasso, mm, by=c("node", "graphNum"))

df3 %>% group_by(measure) %>% summarize(cor(value, poploading))


dd <- demo1()
save(file="psychopathology_networks_demo1_4Nov2016.RData", dd)
str(dd[[1]]$gmetrics$EBICglasso)
options(dplyr.print_max = 1000)
#will be a correlation of nindicators (10) x nreplications (200). So treats indicators as fungible
dd[[3]]$gmetrics$EBICglasso %>% group_by(measure) %>% summarize(cor(value, poploading))
dd[[3]]$gmetrics$EBICglasso %>% group_by(measure) %>% summarize(cor(value, fittedloading))

#aggregate measures (both fitted loadings and centrality) within node before correlation
#then we end up with a correlation between mean loadings and mean centrality measures across replications
dd[[3]]$gmetrics$EBICglasso %>% group_by(measure, node) %>% summarize(mload=mean(fittedloading), mvalue=mean(value)) %>% 
    group_by(measure) %>% summarize(cor(mload,mvalue))

#could imagine the opposite: means of correlations within replications
#so, for each replication, correlate fitted loadings with centrality measure across nodes (10)
#then take the mean of those centrality-loading correlations
dd[[3]]$gmetrics$EBICglasso %>% group_by(graphNum, measure) %>% summarize(cv=cor(fittedloading,value)) %>% 
    group_by(measure) %>% summarize(mcv=mean(cv))

