## For paper: nodal centrality versus factor metrics
## This copies out the sim 1 code from the Rmd file and expands it to more factors.
library(qgraph)
library(igraph)
library(lavaan)
library(viridis)
library(dplyr)
library(lme4)
library(tidyr)
library(broom)
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}

setwd("~/Data_Analysis/psychopathology_network_sim")
source("code/psychopathology_sim_functions.R")

#load("data/onefac_cfa_20examples_200reps_n400.RData")
load("data/onefac_cfa_50examples_100reps_n400.RData")
outstruct <- summarize_loadings_convergence(dd)

png("figures/one-factor EBICGLASSO convergence.png", width=9, height=6, units="in", res=200)
do.call(plot_grid, c(outstruct[["EBICglasso"]]$glist, ncol=1))
dev.off()

png("figures/one-factor pcor convergence.png", width=9, height=6, units="in", res=200)
do.call(plot_grid, c(outstruct[["pcor"]]$glist, ncol=1))
dev.off()

#load("data/twofac_cfa_20examples_200reps_n400.RData")
#load("data/twofac_corr0p4_cfa_50examples_100reps_n400.RData")
load("data/twofac_cfa_50examples_100reps_n400.RData")
outstruct <- summarize_loadings_convergence(dd2)

#mm <- cfa(dd2_corr[[1]]$specification$syntax$fitsyntax, data=dd2_corr[[1]]$simdata[[2]])
mm <- cfa(dd2[[1]]$specification$syntax$fitsyntax, data=dd2[[1]]$simdata[[2]]) #check the model that was fitted
summary(mm,fit.measures=T, standardized=T)
cor(dd2[[1]]$simdata[[2]])

ggplot(filter(dd2[[4]]$gmetrics$EBICglasso, node %in% paste0("y", 1:10)), aes(x=value)) +
  facet_grid(node~measure, scales="free") +
  geom_histogram() + theme_bw(base_size=15)

xx <- filter(dd2[[4]]$gmetrics$EBICglasso, measure=="betweenness" & value > 10)
hist(table(xx$node))

#would be nice to plot the distribution of off-factor associations compared to on-factor

#look at level of off-factor relationship in EBICglasso data
#off_factor <- apply(dd2[[4]]$adjmats$EBICglasso$concat, 1, function(replication) {
#    off_factor <- rep(NA, ncol(replication))
#    for (o in 1:length(off_factor)) {
#      if (o <= 10) {
#        block <- replication[o, 11:20] #this is hard-coded for 10 indicators per factor
#      } else {
#        block <- replication[o, 1:10]
#      }
#      #wvec <- block[lower.tri(block)] #this is wrong because block is a vector
#      off_factor[o] <- sum(abs(block)) #sum(wvec[wvec < 0]) #sum(block[lower.tri(block)])
#    }
#    return(off_factor)
#  })

#library(heatmaply)
##waiting on pandoc
#heatmaply(dd2[[7]]$adjmats$pearson$average, Rowv=FALSE, Colv=FALSE, margins = c(35, 35), column_text_angle = 0, file = "figures/heatmaply_plot.png")
#browseURL("figures/heatmaply_plot.pdf")

#mdf <- reshape2::melt(dd2[[7]]$adjmats$pearson$average)

example_to_plot <- 8
##need a multi-panel plot. two-fac pearson + glasso; two-fac corr pearson + glasso
mdf <- reshape2::melt(dd2[[example_to_plot]]$adjmats$pearson$concat[1,,]) %>% filter(Var1 != Var2)
mdf_glasso <- reshape2::melt(dd2[[example_to_plot]]$adjmats$EBICglasso$concat[1,,]) %>% filter(Var1 != Var2)

#get a graph of the same data
#pdf("figures/two-fac orthogonal glasso.pdf", width=6, height=4)
#qgraph(dd2[[example_to_plot]]$adjmats$EBICglasso$concat[1,,])
#dev.off()

pdf("figures/two-fac orthogonal glasso.pdf", width=10, height=8)
library(ggnetwork)
g_orig <- graph.adjacency(dd2[[example_to_plot]]$adjmats$EBICglasso$concat[1,,], mode="undirected", weighted=TRUE, diag=FALSE)
#n <- ggnetwork(g_orig, layout = "fruchtermanreingold", cell.jitter = 0.75)
n <- ggnetwork(g_orig, layout = "circle", cell.jitter = 0.75)
n <- as.data.frame(lapply(n, as.vector))
ng1 <- ggplot(n, aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_edges(aes(size=weight^3, color=weight)) +
  geom_nodes(color="grey80", size = 21) +
  geom_nodetext(aes(label = vertex.names), size = 8) +
  scale_size(range=c(0.1, 4), guide=FALSE) +
  scale_color_viridis("Partial\ncorrelation", limits=c(-.02, .25), breaks=c(0, .1, .2)) +
  #scale_color_distiller("Partial correlation", palette="RdPu") +
  expand_limits(x=c(-0.05, 1.05), y=c(-0.02, 1.02)) + #avoid the labels going out of bounds
  theme_blank(base_size=16) 
plot(ng1)
dev.off()


pdf("figures/heatmap_pearson_twofactor.pdf", width=8, height=7)
g1 <- ggplot(mdf, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  scale_fill_viridis(expression(paste("Correlation\n(Pearson", italic("r"), ")")), limits=c(-.05, .8), option="inferno") + ylab("") + xlab("Item") + #breaks=c(0, 0.2, 0.4, 0.6)
  #scale_fill_distiller(expression(paste("Correlation\n(Pearson", italic("r"), ")")), limits=c(-.05, 0.8), palette="PuBu") + ylab("") + xlab("Item") + # breaks=c(0.2, 0.4, 0.6)
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5), axis.title.x=element_text(margin=margin(t=15))) + coord_fixed()
plot(g1)
dev.off()

pdf("figures/heatmap_glasso_twofactor.pdf", width=8, height=7)
g2 <- ggplot(mdf_glasso, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  scale_fill_viridis(expression(paste("Correlation\n(Partial", italic("r"), ")")), limits=c(-.02, .25), breaks=c(0, .1, .2)) + ylab("") + xlab("Item") +
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5), axis.title.x=element_text(margin=margin(t=15))) + coord_fixed()
plot(g2)
dev.off()

round(dd2_corr[[7]]$adjmats$pearson$concat[1,,], 2)
round(dd2[[7]]$adjmats$pearson$concat[1,,], 2)

load("data/twofac_corr0p4_cfa_50examples_100reps_n400.RData")

example_to_plot <- 12
mdf <- reshape2::melt(dd2_corr[[example_to_plot]]$adjmats$pearson$concat[1,,]) %>% filter(Var1 != Var2)
mdf_glasso <- reshape2::melt(dd2_corr[[example_to_plot]]$adjmats$EBICglasso$concat[1,,]) %>% filter(Var1 != Var2)

pdf("figures/two-fac correlated glasso.pdf", width=10, height=8)
library(ggnetwork)
g_orig <- graph.adjacency(dd2_corr[[example_to_plot]]$adjmats$EBICglasso$concat[1,,], mode="undirected", weighted=TRUE, diag=FALSE)
#n <- ggnetwork(g_orig, layout = "fruchtermanreingold", cell.jitter = 0.75)
n <- ggnetwork(g_orig, layout = "circle", cell.jitter = 0.75)
n <- as.data.frame(lapply(n, as.vector))
ng2 <- ggplot(n, aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_edges(aes(size=weight^3, color=weight)) +
  geom_nodes(color="grey80", size = 21) +
  geom_nodetext(aes(label = vertex.names), size = 8) +
  scale_size(range=c(0.1, 4), guide=FALSE) +
  scale_color_viridis("Partial\ncorrelation", limits=c(-.02, .25), breaks=c(0, .1, .2)) +
  expand_limits(x=c(-0.05, 1.05), y=c(-0.02, 1.02)) + #avoid the labels going out of bounds
  #scale_color_distiller("Partial correlation", palette="RdPu") +
  theme_blank(base_size=16) 
plot(ng2)
dev.off()

#qgraph(dd2_corr[[example_to_plot]]$adjmats$EBICglasso$concat[1,,])#, legend=TRUE, borders=FALSE, label.cex=0.4, ) #plot average under pure partial correlation
#dev.off()

outstruct <- summarize_loadings_convergence(dd2_corr)

pdf("figures/heatmap_pearson_twofactor_corr.pdf", width=8, height=7)
g3 <- ggplot(mdf, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  scale_fill_viridis(expression(paste("Correlation\n(Pearson", italic("r"), ")")), limits=c(-0.05, 0.8), option="inferno") + ylab("") + xlab("Item") + # breaks=c(0.2, 0.4, 0.6)
  #scale_fill_distiller(expression(paste("Correlation\n(Pearson", italic("r"), ")")), limits=c(-0.05, 0.8), palette="PuBu") + ylab("") + xlab("Item") + # breaks=c(0.2, 0.4, 0.6)
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5), axis.title.x=element_text(margin=margin(t=15))) + coord_fixed()
plot(g3)
dev.off()

pdf("figures/heatmap_glasso_twofactor_corr.pdf", width=8, height=7)
g4 <- ggplot(mdf_glasso, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  scale_fill_viridis(expression(paste("Correlation\n(Partial", italic("r"), ")")), limits=c(-.02, .25), breaks=c(0, .1, .2)) + ylab("") + xlab("Item") +
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5), axis.title.x=element_text(margin=margin(t=15))) + coord_fixed()
plot(g4)
dev.off()

pdf("figures/heatmap_all.pdf", width=15, height=12)
nomargin <- theme(plot.margin=margin(b=5, t=5, l=5, r=5))
plot_grid(g1 + xlab(""), g2 + xlab(""), g3 + xlab(""), g4 + xlab(""), ncol=2, align="hv", labels=c("a)", "b)", "c)", "d)"), label_size=18)
dev.off()

pdf("figures/ebicglasso_networks_both.pdf", width=15, height=6)
nomargin <- theme(plot.margin=margin(b=5, t=5, l=5, r=5))
plot_grid(ng1 + theme(legend.position="none", plot.margin=margin(t=30)), ng2 + theme(plot.margin=margin(t=30)), ncol=2, labels=c("Two-Factor Orthogonal", "Two-Factor Correlated"), 
  label_size=18, label_x=c(0.5, 0.43), hjust=0.5, rel_widths=c(0.85, 1))
dev.off()


load("data/threefac_corr0p4_cfa_50examples_100reps_n400.RData")
mdf <- reshape2::melt(dd3_corr[[7]]$adjmats$pearson$concat[1,,])
mdf_glasso <- reshape2::melt(dd3_corr[[7]]$adjmats$EBICglasso$concat[1,,])

outstruct <- summarize_loadings_convergence(dd3_corr)

load("data/threefac_cfa_50examples_100reps_n400.RData")
mdf <- reshape2::melt(dd3[[7]]$adjmats$pearson$concat[1,,])
mdf_glasso <- reshape2::melt(dd3[[7]]$adjmats$EBICglasso$concat[1,,])

outstruct <- summarize_loadings_convergence(dd3)

#combine all relationships across conditions into one big model

wild_data1 <- do.call(rbind, lapply(1:length(dd), function(e) {
      df <- dd[[e]]$graph_v_factor$EBICglasso$metric_v_loadings
      df$example <- e
      df$nfactors <- 1
      df$corrcond <- "uncorrelated"
      df$method <- "ebicglasso"
      return(df)
    }))

wild_data1_pcor <- do.call(rbind, lapply(1:length(dd), function(e) {
      df <- dd[[e]]$graph_v_factor$pcor$metric_v_loadings
      df$example <- e
      df$nfactors <- 1
      df$corrcond <- "uncorrelated"
      df$method <- "pcor"
      return(df)
    }))

wild_data2 <- do.call(rbind, lapply(1:length(dd2), function(e) {
    df <- dd2[[e]]$graph_v_factor$EBICglasso$metric_v_loadings
    df$example <- e
    df$nfactors <- 2
    df$corrcond <- "uncorrelated"
    df$method <- "ebicglasso"
    return(df)
  }))

wild_data2_pcor <- do.call(rbind, lapply(1:length(dd2), function(e) {
      df <- dd2[[e]]$graph_v_factor$pcor$metric_v_loadings
      df$example <- e
      df$nfactors <- 2
      df$corrcond <- "uncorrelated"
      df$method <- "pcor"
      return(df)
    }))

wild_data2corr <- do.call(rbind, lapply(1:length(dd2_corr), function(e) {
      df <- dd2_corr[[e]]$graph_v_factor$EBICglasso$metric_v_loadings
      df$example <- e
      df$nfactors <- 2
      df$corrcond <- "correlated"
      df$method <- "ebicglasso"
      return(df)
    }))

wild_data2corr_pcor <- do.call(rbind, lapply(1:length(dd2_corr), function(e) {
      df <- dd2_corr[[e]]$graph_v_factor$pcor$metric_v_loadings
      df$example <- e
      df$nfactors <- 2
      df$corrcond <- "correlated"
      df$method <- "pcor"
      return(df)
    }))

wild_data3 <- do.call(rbind, lapply(1:length(dd3), function(e) {
      df <- dd3[[e]]$graph_v_factor$EBICglasso$metric_v_loadings
      df$example <- e
      df$nfactors <- 3
      df$corrcond <- "uncorrelated"
      df$method <- "ebicglasso"
      return(df)
    }))

wild_data3_pcor <- do.call(rbind, lapply(1:length(dd3), function(e) {
      df <- dd3[[e]]$graph_v_factor$pcor$metric_v_loadings
      df$example <- e
      df$nfactors <- 3
      df$corrcond <- "uncorrelated"
      df$method <- "pcor"
      return(df)
    }))

wild_data3corr <- do.call(rbind, lapply(1:length(dd3_corr), function(e) {
      df <- dd3_corr[[e]]$graph_v_factor$EBICglasso$metric_v_loadings
      df$example <- e
      df$nfactors <- 3
      df$corrcond <- "correlated"
      df$method <- "ebicglasso"
      return(df)
    }))

wild_data3corr_pcor <- do.call(rbind, lapply(1:length(dd3_corr), function(e) {
      df <- dd3_corr[[e]]$graph_v_factor$pcor$metric_v_loadings
      df$example <- e
      df$nfactors <- 3
      df$corrcond <- "correlated"
      df$method <- "pcor"
      return(df)
    }))

rm(dd, dd2, dd3, dd2_corr, dd3_corr)

wild_data <- rbind(wild_data1, wild_data1_pcor, wild_data2, wild_data2_pcor, wild_data2corr, wild_data2corr_pcor, wild_data3, wild_data3_pcor, wild_data3corr, wild_data3corr_pcor)
wild_data$nfactors <- factor(wild_data$nfactors)

save(wild_data, file="data/all_loadings_v_graphs_Mar2018.RData")

library(lme4)
library(emmeans)
#these models just get insane and I don't really care about interactions
#go for targeted comparisons below and set aside the NHST mentality
strength <- filter(wild_data, measure=="strength")
mm <- lmer(value ~ fittedloading + nfactors + corrcond + method + (1|example/graphNum) + (1|node), strength)
mm2 <- lmer(value ~ fittedloading + nfactors + corrcond + method + (1|example) + (1|graphNum) + (1|node), strength) #not nested
mm3 <- lmer(value ~ fittedloading + nfactors * method * corrcond + (1|example) + (1|graphNum) + (1|node), strength)
#mm3 <- lmer(value ~ fittedloading * nfactors * corrcond + (1|example/graphNum) + (1|node), strength) #proobably not valid/crucial since there is no correlated one-factor scenario
summary(mm)
emmeans(mm, ~nfactors)
summary(mm2)
anova(mm, mm2)


library(effsize)
#MLMs are super-overpowered here -- will scream about any small change
#simple test 1: is there a decline for strength going from 1 to 2 factors using EBIC and PCOR
load("data/onefac_cfa_50examples_100reps_n400.RData")
outstruct_1 <- summarize_loadings_convergence(dd)
ebic_1 <- outstruct_1$EBICglasso$corrvgraph_detailed %>% ungroup() %>% mutate(nfactors="1")
pcor_1 <- outstruct_1$pcor$corrvgraph_detailed %>% ungroup() %>% mutate(nfactors="1")

load("data/twofac_cfa_50examples_100reps_n400.RData")
outstruct_2 <- summarize_loadings_convergence(dd2)
ebic_2 <- outstruct_2$EBICglasso$corrvgraph_detailed %>% ungroup() %>% mutate(nfactors="2")
pcor_2 <- outstruct_2$pcor$corrvgraph_detailed %>% ungroup() %>% mutate(nfactors="2")

ebic_1v2 <- rbind(ebic_1, ebic_2)
dplyr::filter(ebic_1v2, measure=="strength") %>% group_by(nfactors) %>% summarize(mean(m_overreps), sd(m_overreps)) 
cohen.d(m_overreps ~ nfactors, dplyr::filter(ebic_1v2, measure=="strength"))

dplyr::filter(ebic_1v2, measure=="closeness") %>% group_by(nfactors) %>% summarize(mean(m_overreps), sd(m_overreps)) 
cohen.d(m_overreps ~ nfactors, dplyr::filter(ebic_1v2, measure=="closeness"))

dplyr::filter(ebic_1v2, measure=="betweenness") %>% group_by(nfactors) %>% summarize(mean(m_overreps), sd(m_overreps)) 
cohen.d(m_overreps ~ nfactors, dplyr::filter(ebic_1v2, measure=="betweenness"))

pcor_1v2 <- rbind(pcor_1, pcor_2)
dplyr::filter(pcor_1v2, measure=="strength") %>% group_by(nfactors) %>% summarize(mean(m_overreps), sd(m_overreps)) 
cohen.d(m_overreps ~ nfactors, dplyr::filter(pcor_1v2, measure=="strength"))

dplyr::filter(pcor_1v2, measure=="closeness") %>% group_by(nfactors) %>% summarize(mean(m_overreps), sd(m_overreps)) 
cohen.d(m_overreps ~ nfactors, dplyr::filter(pcor_1v2, measure=="closeness"))

dplyr::filter(pcor_1v2, measure=="betweenness") %>% group_by(nfactors) %>% summarize(mean(m_overreps), sd(m_overreps)) 
cohen.d(m_overreps ~ nfactors, dplyr::filter(pcor_1v2, measure=="betweenness"))


rm(dd, dd2)


s1v2 <- droplevels(dplyr::filter(wild_data, measure=="strength" & method=="ebicglasso" & nfactors %in% c("1", "2")))
#mm1 <- lmer(value ~ fittedloading + nfactors + (1|example/graphNum) + (1|node), s1v2)
cohen.d(value ~ nfactors, s1v2)

s1v2 <- droplevels(dplyr::filter(wild_data, measure=="strength" & method=="pcor" & nfactors %in% c("1", "2")))
cohen.d(value ~ nfactors, s1v2)



s1v2 <- droplevels(dplyr::filter(wild_data, measure=="strength" & method=="pcor" & nfactors %in% c("1", "2")))
mm1 <- lmer(value ~ fittedloading + nfactors + (1|example/graphNum) + (1|node), s1v2)
library(effsize)
cohen.d(value ~ nfactors, s1v2)

cl <- filter(wild_data, measure=="closeness")
mm_cl <- lmer(value ~ fittedloading + nfactors + corrcond + method + (1|example/graphNum) + (1|node), cl)
mm2_cl <- lmer(value ~ fittedloading + nfactors + corrcond + method + (1|example) + (1|graphNum) + (1|node), strength) #not nested

summary(mm)
summary(mm2)
anova(mm, mm2)


#demonstration of off-factor phenomenon
ofd <- get_off_factor_dist(dd2[[4]])
ggplot(filter(ofd$off, value != 0), aes(x=value)) + geom_histogram(bins=30) # + scale_y_continuous(trans='log10')

#what is the ratio of the *summed* off-factor versus *mean* on-factor?
#need to average over replications to avoid divide by 0 problems
mean(ofd$on %>% group_by(node) %>% summarize(value=mean(value)) %>% pull(value) / ofd$off %>% group_by(node) %>% summarize(value=mean(value)) %>% pull(value)) #so mean of on diagonal is 113x larger than off-diagonal

#off_factor <- melt(off_factor, varnames=c("node", "graphNum"), value.name="off_factor_sum") %>% mutate(node=paste0("y", node))
off_factor <- ofd$off %>% mutate(node=paste0("y", node)) %>% rename(off_factor_sum=value)
gmeasures <- dd2[[4]]$gmetrics$EBICglasso %>% inner_join(off_factor, by=c("node", "graphNum"))

#dd2[[4]]$adjmats$EBICglasso$concat[1,,] #spot check

#combine f1, f2 loading columns (since they are only not NA for the primary loading)
gmeasures <- gmeasures %>% gather(key="col", value="fittedloading", f1_fittedloading, f2_fittedloading) %>%
  filter(!is.na(fittedloading))

##Off-factor loadings for ALL EXAMPLES

#need to lapply and combine all examples
get_off_factor_dist_allexamples <- function(exampleslist, edge_method="EBICglasso") {
  ofd_all <- do.call(rbind, lapply(1:length(exampleslist), function(example) {
        
        ofd <- get_off_factor_dist(exampleslist[[example]], method=edge_method)
        ofd$on <- ofd$on %>% dplyr::rename(on_factor_sum=value)
        ofd$off <- ofd$off %>% dplyr::rename(off_factor_sum=value)
        
        ofd_both <- inner_join(ofd$off, ofd$on, by=c("graphNum", "node")) %>% dplyr::mutate(node=paste0("y", node))
        
        gmeasures <- exampleslist[[example]]$gmetrics[[edge_method]] 
        ofd_both <- inner_join(ofd_both, gmeasures, by=c("node", "graphNum"))
        ofd_both$example <- example
        return(ofd_both)
      }))
  
  return(ofd_all)
}

ofd_all_orthogonal_glasso <- get_off_factor_dist_allexamples(dd2, edge_method="EBICglasso")
ofd_all_orthogonal_pcor <- get_off_factor_dist_allexamples(dd2, edge_method="pcor")

#hacky way to get correlated results for two-panel figure (uncomment and re-run, store into ggcorr object...)
ofd_all_orthogonal_glasso <- get_off_factor_dist_allexamples(dd2_corr, edge_method="EBICglasso")
ofd_all_orthogonal_pcor <- get_off_factor_dist_allexamples(dd2_corr, edge_method="pcor")


#combine f1, f2 loading columns (since they are only not NA for the primary loading)
ofg_all_orthogonal_glasso <- ofd_all_orthogonal_glasso %>% gather(key="col", value="fittedloading", f1_fittedloading, f2_fittedloading) %>%
  filter(!is.na(fittedloading)) %>% filter(off_factor_sum > 0) #drop all zero (lassoed) off-factor entries
ofg_all_orthogonal_pcor <- ofd_all_orthogonal_pcor %>% gather(key="col", value="fittedloading", f1_fittedloading, f2_fittedloading) %>%
  filter(!is.na(fittedloading)) %>% filter(off_factor_sum > 0) #for pcor, there shouldn't be any zero entries, but just in case

#so, in short, path metrics are extremely vulnerable to spurious associations with variables that form a subgraph

#divide by metric for MLMs
#betweenness is not log-transformed inside the gmetrics element
st_df <- filter(ofg_all_orthogonal_glasso, measure=="strength") %>% mutate_at(vars(value, fittedloading, off_factor_sum), funs(std=as.vector(scale(.))))
cl_df <- filter(ofg_all_orthogonal_glasso, measure=="closeness") %>% mutate_at(vars(value, fittedloading, off_factor_sum), funs(std=as.vector(scale(.))))
bw_df <- filter(ofg_all_orthogonal_glasso, measure=="betweenness") %>% mutate(value=log10(value + 1.1) + 1) %>% mutate_at(vars(value, fittedloading, off_factor_sum), funs(std=as.vector(scale(.))))
st_df_pcor <- filter(ofg_all_orthogonal_pcor, measure=="strength") %>% mutate_at(vars(value, fittedloading, off_factor_sum), funs(std=as.vector(scale(.))))
cl_df_pcor <- filter(ofg_all_orthogonal_pcor, measure=="closeness") %>% mutate_at(vars(value, fittedloading, off_factor_sum), funs(std=as.vector(scale(.))))
bw_df_pcor <- filter(ofg_all_orthogonal_pcor, measure=="betweenness") %>% mutate(value=log10(value + 1.1) + 1) %>% mutate_at(vars(value, fittedloading, off_factor_sum), funs(std=as.vector(scale(.))))

hist(cl_df$off_factor_sum)
sum(cl_df$off_factor_sum==0)
sum(cl_df$off_factor_sum!=0) #so ~45% have some sort of non-zero effect

hist(cl_df_pcor$off_factor_sum)
sum(cl_df_pcor$off_factor_sum==0)
sum(cl_df_pcor$off_factor_sum!=0) #so ~45% have some sort of non-zero effect

#compute the ratio of off- versus on-factor sums across all datasets
#need to average over replications to avoid divide by 0 problems
mean(st_df %>% group_by(node) %>% summarize(value=mean(on_factor_sum)) %>% pull(value) / st_df %>% group_by(node) %>% summarize(value=mean(off_factor_sum)) %>% pull(value)) #so mean of on diagonal is 113x larger than off-diagonal
sd(st_df %>% group_by(node) %>% summarize(value=mean(on_factor_sum)) %>% pull(value) / st_df %>% group_by(node) %>% summarize(value=mean(off_factor_sum)) %>% pull(value)) #so mean of on diagonal is 113x larger than off-diagonal

head(st_df)


#NB. putting the on_factor_sum as a regressor leads to insane collinearity with fitted loading (not surprising)
#The MLMs should be parameterized in terms of graph and node nested within example.
#This seems right since nodes are largely influenced by population factor loadings, which are drawn at the example level
#Likewise, graphs are largely similar because they are sampled from an example
#Turns out this fits best anyhow

#summary(m_st <- lmer(value ~ fittedloading + off_factor_sum + (1|node) + (1|graphNum) + (1|example), st_df))
#run on (completely) standardized data to get std betas. The function only computes the coefs, not the SEs
#summary(m_st_glasso <- lmer(value ~ fittedloading + off_factor_sum + (1|example/graphNum) + (1|example:node), st_df)) #graphs and nodes are nested within examples
summary(m_st_glasso <- lmer(value_std ~ fittedloading_std + off_factor_sum_std + (1|example/graphNum) + (1|example:node), st_df)) #graphs and nodes are nested within examples
summary(m_st_pcor <- lmer(value_std ~ fittedloading_std + off_factor_sum_std + (1|example/graphNum) + (1|example:node), st_df_pcor)) #graphs and nodes are nested within examples

#confidence intervals, rather than SEs for plotting (otherwise it will just be too tiny to even see)

tmpglasso <- tidy(m_st_glasso) %>% filter(term %in% c("fittedloading_std", "off_factor_sum_std")) %>% mutate(metric="strength", method="ebicglasso") #stdbeta=lm.beta.lmer(m_st_glasso),
cis <- tidy(confint.merMod(m_st_glasso, method="Wald", level=0.99)) %>% filter(.rownames %in% c("fittedloading_std", "off_factor_sum_std")) %>% select(X0.5.., X99.5..) %>%
  rename(ci_0.5=X0.5.., ci_99.5=X99.5..)  #fast computation of CIs given large dataset
tmpglasso <- cbind(tmpglasso, cis)

tmppcor <- tidy(m_st_pcor) %>% filter(term %in% c("fittedloading_std", "off_factor_sum_std")) %>% mutate(metric="strength", method="pcor") #stdbeta=lm.beta.lmer(m_st_pcor),
cis <- tidy(confint.merMod(m_st_pcor, method="Wald", level=0.99)) %>% filter(.rownames %in% c("fittedloading_std", "off_factor_sum_std")) %>% select(X0.5.., X99.5..) %>%
  rename(ci_0.5=X0.5.., ci_99.5=X99.5..)  #fast computation of CIs given large dataset
tmppcor <- cbind(tmppcor, cis)

par_st <- rbind(tmpglasso, tmppcor)
#consider only case where there is a non-zero sum
#summary(m_st <- lmer(value ~ fittedloading + off_factor_sum + (1|example/graphNum) + (1|example:node), filter(st_df, off_factor_sum > 0))) #graphs and nodes are nested within examples

summary(m_cl_glasso <- lmer(value_std ~ fittedloading_std + off_factor_sum_std + (1|example/graphNum) + (1|example:node), cl_df)) #graphs and nodes are nested within examples
summary(m_cl_pcor <- lmer(value_std ~ fittedloading_std + off_factor_sum_std + (1|example/graphNum) + (1|example:node), cl_df_pcor)) #graphs and nodes are nested within examples

tmpglasso <- tidy(m_cl_glasso) %>% filter(term %in% c("fittedloading_std", "off_factor_sum_std")) %>% mutate(metric="closeness", method="ebicglasso") #stdbeta=lm.beta.lmer(m_st_glasso),
cis <- tidy(confint.merMod(m_cl_glasso, method="Wald", level=0.99)) %>% filter(.rownames %in% c("fittedloading_std", "off_factor_sum_std")) %>% select(X0.5.., X99.5..) %>%
  rename(ci_0.5=X0.5.., ci_99.5=X99.5..)  #fast computation of CIs given large dataset
tmpglasso <- cbind(tmpglasso, cis)

tmppcor <- tidy(m_cl_pcor) %>% filter(term %in% c("fittedloading_std", "off_factor_sum_std")) %>% mutate(metric="closeness", method="pcor") #stdbeta=lm.beta.lmer(m_st_pcor),
cis <- tidy(confint.merMod(m_cl_pcor, method="Wald", level=0.99)) %>% filter(.rownames %in% c("fittedloading_std", "off_factor_sum_std")) %>% select(X0.5.., X99.5..) %>%
  rename(ci_0.5=X0.5.., ci_99.5=X99.5..)  #fast computation of CIs given large dataset
tmppcor <- cbind(tmppcor, cis)

par_cl <- rbind(tmpglasso, tmppcor)

summary(m_bw_glasso <- lmer(value_std ~ fittedloading_std + off_factor_sum_std + (1|example/graphNum) + (1|example:node), bw_df)) #graphs and nodes are nested within examples
summary(m_bw_pcor <- lmer(value_std ~ fittedloading_std + off_factor_sum_std + (1|example/graphNum) + (1|example:node), bw_df_pcor)) #graphs and nodes are nested within examples

tmpglasso <- tidy(m_bw_glasso) %>% filter(term %in% c("fittedloading_std", "off_factor_sum_std")) %>% mutate(metric="betweenness", method="ebicglasso") #stdbeta=lm.beta.lmer(m_st_glasso),
cis <- tidy(confint.merMod(m_bw_glasso, method="Wald", level=0.99)) %>% filter(.rownames %in% c("fittedloading_std", "off_factor_sum_std")) %>% select(X0.5.., X99.5..) %>%
  rename(ci_0.5=X0.5.., ci_99.5=X99.5..)  #fast computation of CIs given large dataset
tmpglasso <- cbind(tmpglasso, cis)

tmppcor <- tidy(m_bw_pcor) %>% filter(term %in% c("fittedloading_std", "off_factor_sum_std")) %>% mutate(metric="betweenness", method="pcor") #stdbeta=lm.beta.lmer(m_st_pcor),
cis <- tidy(confint.merMod(m_bw_pcor, method="Wald", level=0.99)) %>% filter(.rownames %in% c("fittedloading_std", "off_factor_sum_std")) %>% select(X0.5.., X99.5..) %>%
  rename(ci_0.5=X0.5.., ci_99.5=X99.5..)  #fast computation of CIs given large dataset
tmppcor <- cbind(tmppcor, cis)

par_bw <- rbind(tmpglasso, tmppcor)

par_all <- rbind(par_st, par_cl, par_bw) %>% mutate(term=recode(term, fittedloading_std="loading", off_factor_sum_std="off-factor\nsum"), 
  metric=ordered(metric, levels=c("strength", "closeness", "betweenness"), labels=c("Strength", "Closeness", "Betweenness")),
  method=if_else(method=="ebicglasso", "EBICGLASSO", "PCor"))
  
  
pdf("figures/loadings_v_off_factor_lmercoefs_orthogonal.pdf", width=8, height=5)
ggorthogonal <- ggplot(par_all, aes(x=term, y=estimate, ymin=ci_0.5, ymax=ci_99.5, color=method, fill=method)) + 
  geom_col(position=position_dodge(width=0.5), width=0.4) + geom_linerange(position=position_dodge(width=0.5), size=2, color="black") + facet_grid(. ~ metric, scales="free") + theme_bw(base_size=14) +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + ylab("Standardized coefficient") + xlab("Predictor") +
  #scale_fill_brewer("Edge definition", palette="Set2") + scale_color_brewer("Edge definition", palette="Set2")
  scale_fill_viridis("Edge definition", discrete=TRUE, begin=0.2, end=0.8, option="cividis") + scale_color_viridis("Edge definition", discrete=TRUE, begin=0.2, end=0.8, option="cividis")

plot(ggorthogonal)
dev.off()

pdf("figures/loadings_v_off_factor_lmercoefs_both.pdf", width=8, height=8)
plot_grid(ggorthogonal, ggcorr + theme(plot.margin=margin(t=30)), nrow=2, align="hv", labels=c("a)", "b)"))
dev.off()

#some leftovers from model fitting
#summary(m_cl <- lmer(value ~ fittedloading + off_factor_sum + (1|node) + (1|graphNum) + (1|example), cl_df))
#summary(m_cl2 <- lmer(value ~ fittedloading + off_factor_sum + (1|node) + (1|example/graphNum), cl_df)) #graphs are nested within examples
#summary(m_cl3 <- lmer(value ~ fittedloading + off_factor_sum + (1|example/graphNum), cl_df)) #graphs are nested within examples
#summary(m_cl4 <- lmer(value ~ fittedloading + off_factor_sum + (1|example/graphNum) + (1|example:node), cl_df)) #graphs and nodes are nested within examples
#summary(m_cl4 <- lmer(value ~ fittedloading + off_factor_sum + (1|example/graphNum) + (1|example:node), filter(cl_df, off_factor_sum > 0))) #graphs and nodes are nested within examples
#summary(m_cl4 <- lmer(value ~ fittedloading + off_factor_sum + (1|example/graphNum) + (1|example:node), filter(cl_df_pcor, off_factor_sum > 0))) #graphs and nodes are nested within examples
#anova(m_cl, m_cl2, m_cl3, m_cl4) #also turns out it clearly fits best
#cor(select(st_df, fittedloading, off_factor_sum, on_factor_sum))
#lm.beta.lmer(m_cl4)
#lm.beta.lmer(m_cl2)
#lm.beta.lmer(m_cl)

#non MLM variant (flawed)
#summary(m_cl2 <- lm(value ~ fittedloading + off_factor_sum + node + factor(example) + factor(graphNum), cl_df)) #graphs are nested within examples

#betweenness
#summary(m_bw1 <- lmer(value ~ fittedloading + off_factor_sum + (1|node) + (1|graphNum) + (1|example), filter(ofd_all, measure=="betweenness" & off_factor_sum > 0)))
#summary(m_bw2 <- lmer(value ~ fittedloading + off_factor_sum + I(off_factor_sum^2) + (1|node) + (1|graphNum), filter(gmeasures, measure=="betweenness" & off_factor_sum > 0)))
#summary(m_bw3 <- lmer(value ~ fittedloading + sqrt(off_factor_sum) + (1|node) + (1|graphNum), filter(gmeasures, measure=="betweenness" & off_factor_sum > 0)))
#summary(m_bw4 <- lmer(value ~ fittedloading + I(1/off_factor_sum) + (1|node) + (1|graphNum), filter(gmeasures, measure=="betweenness" & off_factor_sum > 0)))



#NOT USING ANYTHING BELOW HERE FOR SIM1 IN PAPER
anova(m_bw1, m_bw2, m_bw3, m_bw4)

hist(sqrt(filter(gmeasures, measure=="betweenness" & off_factor_sum > 0)$off_factor_sum))

##conclusion: closeness and betweenness are *strongly* influenced by the chance covariation of indicators with the other factor
ggplot(filter(gmeasures, measure=="betweenness" & off_factor_sum > 0), aes(x=sqrt(off_factor_sum), y=value, color=fittedloading)) +
  geom_point() + stat_smooth(method="lm")

filter(gmeasures, measure=="betweenness" & off_factor_sum > 0) %>% summarize(cor(off_factor_sum, value))

summary(lm(log10(value +1.1) ~ fittedloading + off_factor_sum, filter(gmeasures, measure=="betweenness" & off_factor_sum > 0)))
summary(lm(value ~ fittedloading + off_factor_sum, filter(gmeasures, measure=="closeness" & off_factor_sum > 0)))
summary(lm(value ~ fittedloading + off_factor_sum, filter(gmeasures, measure=="strength")))
summary(lm(value ~ fittedloading, filter(gmeasures, measure=="strength")))

#OLDER STUFF

png("figures/two-factor EBICGLASSO convergence.png", width=9, height=6, units="in", res=200)
do.call(plot_grid, c(outstruct[["EBICglasso"]]$glist, ncol=1))
dev.off()

png("figures/two-factor pcor convergence.png", width=9, height=6, units="in", res=200)
do.call(plot_grid, c(outstruct[["pcor"]]$glist, ncol=1))
dev.off()



ggplot(filter(dd[[1]]$graph_v_factor$EBICglasso$metric_v_loadings, graphNum==1), aes(x=fittedloading, y=value)) +
  geom_point() + facet_wrap(~measure, scales="free")


kable(corrvgraph_pcor_detailed %>% group_by(measure, factor) %>% summarize_all(funs(mean, sd)), 
  digits=4, caption = "Association of PCOR graph measures with 1-factor CFA loadings")

smat <- do.call(rbind, lapply(dd, function(rep) {
      brep <- filter(rep$graph_v_factor$EBICglasso$metric_v_loadings, node=="y1" & measure=="strength" & fittedloading < 1)
    }))

#smat$graphNum <- rep(1:(nrow(smat)/2), each=2) #smat has two rows per replication. Need to identify these to get the spread to work properly

bb <- smat %>% select(graphNum, value, factor, fittedloading) %>% spread(key=factor, value=fittedloading) #graphNum
g <- ggplot(bb, aes(x=f1, y=value)) + geom_point(alpha=0.5) + stat_smooth(method="lm") + xlab("Fitted factor loading") +
  ylab("Strength centrality") + annotate(geom="text", x = 0.4, y=1.1, label="r = 0.97") + theme_cowplot(font_size=20)

pdf("figures/strength 1-factor loadings.pdf", width=5, height=4)
plot(g)
dev.off()


kable(
  broom::tidy(cor.test(~ fittedloading + value, filter(smat, factor=="f1"))), 
  digits=2, caption="Association of strength with primary factor loading (f1)")


cmat <- do.call(rbind, lapply(dd, function(rep) {
      brep <- filter(rep$graph_v_factor$EBICglasso$metric_v_loadings, node=="y1" & measure=="closeness" &
          fittedloading < 1)
    }))

cmat$graphNum <- rep(1:(nrow(cmat)/2), each=2) #smat has two rows per replication. Need to identify these to get the spread to work properly

bb <- cmat %>% select(graphNum, value, factor, fittedloading) %>% spread(key=factor, value=fittedloading) #graphNum
g <- ggplot(bb, aes(x=f1, y=value)) + geom_point(alpha=0.3) + stat_smooth(method="lm") + xlab("Fitted factor loading") + ylab("Closeness centrality") + annotate(geom="text", x = 0.4, y=0.15, label="r = 0.96") +
  theme_cowplot(font_size=20)
pdf("figures/closeness 1-factor loadings.pdf", width=5, height=4)
plot(g)
dev.off()