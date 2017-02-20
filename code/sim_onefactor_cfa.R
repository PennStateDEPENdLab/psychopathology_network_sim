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
library(igraph)

x <- demomaster(nexamples=1, nindicators=24, nfactors=3, nreplications=200, n=400, 
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

x <- demomaster(nexamples=1, nindicators=24, nfactors=1, nreplications=200, n=400, 
    errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

x <- demomaster(nexamples=1, nindicators=24, nfactors=1, nreplications=200, n=400, 
    errorvars="eqvars", ivar=1, thetacorlist=NULL)

x[[1]]$graph_v_factor$EBICglasso$corr_v_fitted

#does transforming betweenness linearize/increase the correlation? Scatterplots suggest a nonlinear relationship
#yes... and this is supported in some literature: https://books.google.com/books?id=rbxPm93PRY8C&pg=PA175&lpg=PA175&dq=transform+betweenness+log&source=bl&ots=Mxly7a-oJ4&sig=Ulo6OCupLaWJU25gT6Y7jbPCsdY&hl=en&sa=X&ved=0ahUKEwj8i6nKjurQAhUN92MKHegwAyEQ6AEIHDAA#v=onepage&q=transform%20betweenness%20log&f=false
with(filter(x[[1]]$graph_v_factor$EBICglasso$metric_v_loadings, measure=="betweenness"),
    cor(fittedloading, log10(value + 1.1)*100))

with(filter(x[[1]]$graph_v_factor$EBICglasso$metric_v_loadings, measure=="betweenness"),
    cor(fittedloading, value, method="spearman"))

filter(x[[1]]$graph_v_factor$EBICglasso$metric_v_loadings, measure=="betweenness") %>%
    ggplot(aes(x=log10(value+1.1)*100, y=fittedloading)) + geom_point() + stat_smooth() #+ scale_x_log2()

#MNH Nov2016: checked whether demomaster supersedes demo1. yes, confirmed that all data structures and results are identical
##First demonstration: correlation of one-factor CFA loadings with centrality measures

dd <- demomaster(nexamples=2, nindicators=10, nfactors=1, nreplications=200, n=400, 
    loadingsampler=seq(.4, .95, .05), errorvars="eqresidvar", ivar=0.5, thetacorlist=NULL)

dd[[1]]$graph_v_factor$EBICglasso$corr_v_fitted
dd[[2]]$graph_v_factor$EBICglasso$corr_v_fitted

dd[[1]]$graph_v_factor$pcor$corr_v_fitted
dd[[2]]$graph_v_factor$pcor$corr_v_fitted



#save(file="psychopathology_networks_demo1_4Nov2016.RData", dd)
#load(file="psychopathology_networks_demo1_4Nov2016.RData")

dd[[1]]$adjmats$EBICglasso$average
dd[[1]]$graph_v_factor$pcor$corr_v_fitted
dd[[1]]$graph_v_factor$EBICglasso$corr_v_fitted
#ebicavg <- list(graph=dd[[1]]$adjmats$EBICglasso$average)
#class(ebicavg) <- c("bootnetResult", "list")

g <- igraph::graph_from_adjacency_matrix(dd[[1]]$adjmats[["EBICglasso"]]$concat[33,,], mode="undirected", weighted=TRUE, diag=FALSE)
g <- delete.edges(g, which(E(g)$weight < 0.01))
cluster_louvain(g)
cluster_fast_greedy(g)
cluster_leading_eigen(g)


pdf("figures/est ebic avg network.pdf", width=8, height=8)
#just plot the average adjacency directly
qgraph(dd[[1]]$adjmats$EBICglasso$average)
dev.off()

#fake <- list(graph=ebicavg) #doesn't work
#plot(fake, layout="spring", labels=TRUE)

options(dplyr.print_max = 1000)
#will be a correlation of nindicators (10) x nreplications (200). So treats indicators as fungible
dd[[2]]$gmetrics$EBICglasso %>% group_by(measure) %>% summarize(cor(value, poploading))
dd[[2]]$gmetrics$EBICglasso %>% group_by(measure) %>% summarize(cor(value, fittedloading))

dd[[2]]$gmetrics$pcor %>% group_by(measure) %>% summarize(cor(value, poploading))
dd[[2]]$gmetrics$pcor %>% group_by(measure) %>% summarize(cor(value, fittedloading))

#aggregate measures (both fitted loadings and centrality) within node before correlation
#then we end up with a correlation between mean loadings and mean centrality measures across replications
dd[[2]]$gmetrics$EBICglasso %>% group_by(measure, node) %>% summarize(mload=mean(fittedloading), mvalue=mean(value)) %>% 
    group_by(measure) %>% summarize(cor(mload,mvalue))

#could imagine the opposite: means of correlations within replications
#so, for each replication, correlate fitted loadings with centrality measure across nodes (10)
#then take the mean of those centrality-loading correlations
dd[[2]]$gmetrics$EBICglasso %>% group_by(graphNum, measure) %>% summarize(cv=cor(fittedloading,value)) %>% 
    group_by(measure) %>% summarize(mcv=mean(cv))

pdf("figures/graph versus loading scatter.pdf", width=7, height=10)
dd[[1]]$gmetrics$EBICglasso %>% group_by(graphNum, measure) %>% 
    ggplot(aes(x=fittedloading, y=value)) + facet_wrap(~measure, ncol=1) + geom_point(alpha=0.5) + theme_bw(base_size=20)
dev.off()

#try to plot the centrality measures and fitted loadings side by side (essentially extending the default centralityPlot)
str(dd[[2]]$gmetrics$EBICglasso)

#need to spread out centrality measures, then join back up
xx <- dd[[2]]$gmetrics$EBICglasso %>% spread(key=measure, value=value) %>% select(-graph, -type, -poploading) %>%
    gather(key=key, value=value, -node, -graphNum) %>% group_by(key) %>% mutate(zval=as.vector(scale(value))) %>% ungroup()


xx$node <- ordered(xx$node, levels=paste0("y", 1:10))
xx$key <- ordered(xx$key, levels=c("fittedloading", "Strength", "Closeness", "Betweenness"), labels=c("Factor Loading", "Strength", "Closeness", "Betweenness"))

#ggplot(select(xx, -poploading), aes(x=node, y=value)) + facet_wrap(~key) + geom_boxplot()
ggplot(xx, aes(x=node, y=value)) + facet_grid(key~., scales="free") + geom_boxplot() + coord_flip()

#ggplot(xx, aes(x=node, y=value)) + facet_grid(key~., scales="free") + stat_summary(fun.data = mean_cl_boot) + coord_flip()
ggplot(xx, aes(x=node, y=zval)) + facet_grid(key~., scales="free") + stat_summary(fun.data = mean_cl_boot) + coord_flip()

pdf("figures/demo1_candidate_figure.pdf", width=7, height=11)
ggplot(xx, aes(x=key, y=zval)) + facet_grid(node~.) + stat_summary(fun.data = mean_cl_boot) + coord_flip() +
    ggtitle("200 replications (n=400)") + theme_bw(base_size=18) + 
    scale_x_discrete(limits = rev(levels(xx$key))) +
    xlab("") + ylab("Standardized measure (z +/- 95% CI)")
dev.off()

xx %>% group_by(key) %>% summarize(m=mean(value), low=min(value), high=max(value))

##DEMO ONE-FACTOR PART 2

##okay, what if we introduce residual correlation but keep everything else the same?
thetacorlist <- list(
    list("y1", "y2", 0.5),
    list("y1", "y3", 0.4),
    list("y4", "y6", 0.5),
    list("y7", "y9", 0.3)
)
dd2 <- demo1(thetacorlist=thetacorlist)

pdf("figures/est ebic avg network error cor.pdf", width=8, height=8)
#just plot the average adjacency directly
qgraph(dd2[[1]]$adjmats$EBICglasso$average)
dev.off()

dd2[[1]]$gmetrics$EBICglasso %>% group_by(measure, node) %>% summarize(mload=mean(fittedloading), mvalue=mean(value)) %>% 
    group_by(measure) %>% summarize(cor(mload,mvalue))

dd2[[1]]$gmetrics$EBICglasso %>% group_by(graphNum, measure) %>% summarize(cv=cor(fittedloading,value)) %>% 
    group_by(measure) %>% summarize(mcv=mean(cv))

#scatter plot of the same
dd2[[1]]$gmetrics$EBICglasso %>% group_by(graphNum, measure) %>% 
    ggplot(aes(x=fittedloading, y=value)) + facet_wrap(~measure) + geom_point()

xx <- dd[[2]]$gmetrics$EBICglasso %>% spread(key=measure, value=value) %>% select(-graph, -type) %>%
    gather(key=key, value=value, -node, -graphNum) %>% group_by(key) %>% mutate(zval=as.vector(scale(value))) %>% ungroup()

xx$node <- ordered(xx$node, levels=paste0("y", 1:10))
xx$key <- ordered(xx$key, levels=c("fittedloading", "poploading", "Strength", "Closeness", "Betweenness"), labels=c("Factor Loading", "Pop. Loading", "Strength", "Closeness", "Betweenness"))
ggplot(xx, aes(x=node, y=zval)) + facet_grid(key~., scales="free") + stat_summary(fun.data = mean_cl_boot) + coord_flip()

#pdf("figures/demo1_candidate_figure_errocor.pdf", width=7, height=11)
#ggplot(xx, aes(x=key, y=zval)) + facet_grid(node~.) + stat_summary(fun.data = mean_cl_boot) + coord_flip() +
#    ggtitle("200 replications (n=400)") + theme_bw(base_size=18) + 
#    scale_x_discrete(limits = rev(levels(xx$key))) +
#    xlab("") + ylab("Standardized measure (z +/- 95% CI)")
#dev.off()

