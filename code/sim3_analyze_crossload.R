library(qgraph)
library(igraph)
library(lavaan)
library(viridis)
library(dplyr)
library(lme4)
library(ggplot2)
library(cowplot)
library(caTools)
library(tidyr)
library(lm.beta)

setwd(file.path(getMainDir(), "psychopathology_network_sim"))
source("code/psychopathology_sim_functions.R")

load("data/sim3_crossload_mlist.RData")

#check the inversion of weights idea against bootnet to make sure I'm not departing from the typical netwerker approach
#vv <- bootnet_EBICglasso(mlist_crossload[[1]]$simdata[[1]])
#
#gg <- igraph::graph_from_adjacency_matrix(vv$graph, mode="undirected", weighted=TRUE, diag=FALSE)
#E(gg)$weight <- abs(E(gg)$weight) #this is how qgraph solves negative edges (adopt for now)
#
#centralityTable(vv$graph, standardized=F) %>% filter(measure=="Strength")
#strength(gg)
#
#centralityTable(vv$graph, standardized=F) %>% filter(measure=="Betweenness")
#betweenness(gg, weights=1/E(gg)$weight) #yes, matches bootnet
#
#centralityTable(vv$graph, standardized=F) %>% filter(measure=="Closeness")
#closeness(gg, weights=1/E(gg)$weight) #yes, matches bootnet
#
#all looks good

#this is the underlying simulation grid for the models
loading_grid <- expand.grid(f1=seq(0.2, 0.8, .05), f2=seq(0.2, 0.8, 0.05)) %>% mutate(condition=1:n())

#curious about 

pdf("figures/sim3 loading grid.pdf", width=7, height=5)
ggplot(loading_grid, aes(x=f1, y=f2, fill=f1 + f2)) + geom_tile() +
  #scale_color_gradient() + 
  scale_fill_viridis("Summed\nloadings", option="viridis", breaks=c(0.4, 0.8, 1.2, 1.6)) +
  xlab("Factor 1 loading") + ylab("Factor 2 loading") +
  theme_cowplot(font_size=20)
dev.off()

#curious about EFA being able to "get it"
df <- mlist_crossload[[117]]$simdata[[1]]
library(psych)
ef1 <- fa(df, fm="pa", nfactors=2, rotate="oblimin")
ef2 <- fa(df, fm="pa", nfactors=2, rotate="oblimin")
vss(df, fm="pa", rotate="oblimin") #yes, gets 2 fac

loading_grid[117,]
cf1 <- 'f1 =~ NA*y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y17
f2 =~ NA*y9 + y10 + y11 + y12 + y13 + y14 + y15 + y16 + y17
f1 ~~ 1.0*f1
f2 ~~ 1.0*f2
f1 ~~ f2
'
mm <- cfa(cf1, df)
summary(mm, standardized=TRUE, fit.measures=TRUE)
modificationIndices(mm, min=3) #yes, clearly sees the missing factor


allmetrics_rand <- do.call(rbind, lapply(1:length(mlist_crossload), function(i) {
        df <- filter(mlist_crossload[[i]]$graph_v_factor$EBICglasso$metric_v_loadings, node=="y17")
        df$condition <- i
        return(df)
      })) %>% arrange(poploading, measure, factor, graphNum) %>%
  arrange(poploading, measure, factor, graphNum) %>%
  select(condition, node, measure, graphNum, value, factor, poploading) %>% spread(key=factor, value=poploading) %>%
  mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)

#add a comparator node, y3, which loads on f1, and y10, which loads on f2
#trying to get multiple DVs (fitted loading and graph metric) is a pain. Go to older melt + dcast
#allmetrics_rand_comparator <- do.call(rbind, lapply(1:length(mlist_crossload), function(i) {
#      df <- filter(mlist_crossload[[i]]$graph_v_factor$EBICglasso$metric_v_loadings, node %in% c("y3", "y10"))
#      df$condition <- i
#      return(df)
#    })) %>% unite(col="node_measure", node, measure) %>% select(-factor, -poploading, -fittedloading) %>% spread(key=node_measure, value=value) #%>% mutate(id=1:n()) 

allmetrics_rand_comparator <- do.call(rbind, lapply(1:length(mlist_crossload), function(i) {
        df <- filter(mlist_crossload[[i]]$graph_v_factor$EBICglasso$metric_v_loadings, node %in% c("y3", "y10"))
        df$condition <- i
        return(df)
      })) %>% select(-factor, -poploading) %>% dplyr::rename(gmetric=value) %>% 
  reshape2::melt(id.vars=c("node", "measure", "graphNum", "condition")) %>% 
  reshape2::dcast(graphNum + condition + measure ~ node + variable)

allmetrics_rand <- allmetrics_rand %>% inner_join(allmetrics_rand_comparator, by=c("graphNum", "measure", "condition"))
smat_wide <- allmetrics_rand %>% filter(measure=="strength") 
cmat_wide <- allmetrics_rand %>% filter(measure=="closeness")
bmat_wide <- allmetrics_rand %>% filter(measure=="betweenness")

#cmat_rand$graphNum <- rep(1:(nrow(cmat_rand)/2), each=2) #smat has two rows per replication. Need to identify these to get the spread to work properly

#look at how the correlation matrix is affected
corr_rand <- do.call(rbind, lapply(mlist_crossload, function(rep) {
      #old code: this approach replicates the average association from the *average* matrix, rather than keeping all replication variability
      #but the f1loading, f2loading below was keeping each replication individually, leading to an odd combined matrix (some variables repeated for all replications, others varying)
      #instead, use the concat field to keep the replication variability
      #wif1 <- rep$adjmats$pearson$average[c(1,3:9), c(1,3:9)]
      #wif1m <- mean(wif1[lower.tri(wif1)])
      #wif2 <- rep$adjmats$pearson$average[10:18, 10:18]
      #wif2m <- mean(wif2[lower.tri(wif2)])
      #bwf1f2m <- mean(rep$adjmats$pearson$average[c(1,3:9), 10:18])
      #f1y17m <- mean(rep$adjmats$pearson$average[2, c(1,3:9)]) #correlation of target item (y17) with other items on f1
      #f2y17m <- mean(rep$adjmats$pearson$average[2, 10:18]) #correlation of target item (y17) with other items on f2
      
      wif1m <- apply(rep$adjmats$pearson$concat, 1, function(x) { subx <- x[c(1,3:9), c(1,3:9)]; mean(subx[lower.tri(subx)]) }) #correlation among items loading on f1
      wif2m <- apply(rep$adjmats$pearson$concat, 1, function(x) { subx <- x[10:18, 10:18]; mean(subx[lower.tri(subx)]) }) #correlation among items loading on f2
      bwf1f2m <- apply(rep$adjmats$pearson$concat, 1, function(x) { subx <- x[c(1,3:9), 10:18]; mean(subx[lower.tri(subx)]) }) #correlation between items on f1 and f2 (ignoring the target)
      f1y17m <- apply(rep$adjmats$pearson$concat, 1, function(x) { subx <- x[2, c(1,3:9)]; mean(subx[lower.tri(subx)]) }) #correlation of target item (y17) with other items on f1
      f2y17m <- apply(rep$adjmats$pearson$concat, 1, function(x) { subx <- x[2, 10:18]; mean(subx[lower.tri(subx)]) }) #correlation of target item (y17) with other items on f2
      
      f1loading <- dplyr::filter(rep$graph_v_factor$EBICglasso$metric_v_loadings, factor=="f1" & node=="y17" & measure=="betweenness") %>% arrange(graphNum) %>% pull(fittedloading) #just filtering to betweenness because loadings are repeated for each metric. yields nrow = nreps (100)
      f2loading <- dplyr::filter(rep$graph_v_factor$EBICglasso$metric_v_loadings, factor=="f2" & node=="y17" & measure=="betweenness") %>% arrange(graphNum) %>% pull(fittedloading) #just filtering to betweenness because loadings are repeated for each metric. yields nrow = nreps (100)
      data.frame(graphNum=1:length(f1loading), wif1m=wif1m, wif2m=wif2m, bwf1f2m=bwf1f2m, f1y17m=f1y17m, f2y17m=f2y17m, f1loading=f1loading, f2loading=f2loading)
#        f1loading=select(rep[["simsemout"]]@coef, matches("f1=~y17")), #no longer saved in out struct
#        f2loading=select(rep[["simsemout"]]@coef, matches("f2=~y17")))
    }))

Hmisc::rcorr(as.matrix(select(corr_rand, -graphNum)))

library(apaTables)
#this is pretty handy as a start
apa.cor.table(as.matrix(select(corr_rand, -graphNum)), filename="/Users/mnh5174/TresorSync/Manuscripts/Hallquist Cross-Sectional Networks Critique 2016/sim3 cross-load means and corrs.doc")

#Conclusions about effects of cross-loading on correlations among items:
#  
#1) Higher factor loading of y17 on f1 goes with higher average correlation (r = .91) of y17 with y1-y9 (other f1 items): obvious
#2) Converse: Higher factor loading of y17 on f2 goes with higher average correlation (r = .91) of y17 with y10-y18 (other f2 items): obvious
#3) Higher factor loading of y17 on f2 goes with lower average correlation (r = -.35) with f1 items, and vice versa (r = -.34). So loading drives similarity to primary factor and dissimilarity to cross-loaded factor.
#4) Non-trivial: Higher correlation of y17 with f1 items associated with lower correlation with f2 items, r = -.63.

sapply(corr_rand, mean)

##STRENGTH
#smat_wide <- smat_rand %>% select(graphNum, value, factor, fittedloading) %>% spread(key=factor, value=fittedloading) #use fitted versus population loading
#smat_wide <- smat_rand %>% select(condition, graphNum, value, factor, poploading) %>% spread(key=factor, value=poploading) %>% #graphNum contains the identifying column
#  mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)
#ggplot(smat_wide, aes(x=f1, y=value)) + geom_point() + stat_smooth()
summary(lm(value ~ f1c + f2c + f1sq + f2sq, smat_wide)) #so measure becomes hugely predicted by each loading
summary(lmer(value ~ f1c + f2c + f1sq + f2sq + (1|condition) + (1|graphNum), smat_wide)) #so measure becomes hugely predicted by each loading

#aggregate replications within each combination of f1 and f2loading
#strength_agg <- smat_wide %>% group_by(f1, f2) %>% summarize(value=mean(value)) %>% ungroup() %>% mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)
strength_agg <- smat_wide %>% group_by(f1, f2) %>% summarize_at(vars(value, y10_gmetric, y3_gmetric, y10_fittedloading, y3_fittedloading), funs(mean)) %>% ungroup() %>% mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)
summary(m1 <- lm(value ~ f1c + f2c, strength_agg)) #so measure becomes hugely predicted by each loading
summary(lm.beta(m1))
summary(m2 <- lm(value ~ f1c + f2c + f1sq + f2sq, strength_agg)) #so measure becomes hugely predicted by each loading
anova(m1, m2)

summary(m1 <- lm(y10_gmetric ~ y10_fittedloading , strength_agg))
summary(lm(y10_gmetric ~ f1c * f2c, strength_agg)) #for a comparator, still driving the bus...
summary(m2 <- lm(y10_gmetric ~ f1c * f2c + y10_fittedloading, strength_agg)) #for a comparator, still driving the bus...
anova(m1, m2)

summary(lm(y10_gmetric ~ f1c * f2c + f1sq + f2sq, strength_agg)) #for a comparator, still driving the bus...
summary(m3 <- lm(y10_gmetric ~ f1c * f2c + f1sq + f2sq + y10_fittedloading , strength_agg))

anova(m1, m2, m3)
with(strength_agg, cor.test(y10_gmetric, value))


#yes, this is interesting....
ggplot(smat_wide, aes(x=f1, y=value, color=f2, group=f2)) + stat_smooth()

pdf("figures/strength f1f2 2018.pdf", width=7, height=5)
g1 <- ggplot(strength_agg, aes(x=f1, y=f2, fill=value)) + geom_tile() +
  #scale_color_gradient() + 
  scale_fill_viridis("Strength", option="viridis", breaks=c(0.5, 1, 1.5)) +
  xlab("y17 factor 1 loading") + ylab("y17 factor 2 loading") +
  theme_cowplot(font_size=20) + coord_fixed()
#ggtitle("Association of f1 (x) and f2 (y) loadings with strength (color)")
plot(g1)
dev.off()

#so, the residue in average strength is largely driven by the condition
pdf("figures/strength f1f2 2018 y10 comparator.pdf", width=7, height=5)
g2 <- ggplot(strength_agg, aes(x=f1, y=f2, fill=y10_gmetric)) + geom_tile() +
  #scale_color_gradient() + 
  scale_fill_viridis("Strength", option="viridis", breaks=c(0.94, 0.98, 1.02)) +
  xlab("y17 factor 1 loading") + ylab("y17 factor 2 loading") +
  theme_cowplot(font_size=20) + coord_fixed()
plot(g2)
#ggtitle("Association of f1 (x) and f2 (y) loadings with strength (color)")
dev.off()


#look at effects on other items
#fortunately,  the fitted loadings are not contaminated by the condition (this is not strength-specific)
summary(lm(y10_fittedloading ~ f1c * f2c + f1sq + f2sq, smat_wide))
summary(lm(y10_fittedloading ~ f1c * f2c + f1sq + f2sq, smat_wide))
summary(lm(y3_fittedloading ~ f1c * f2c + f1sq + f2sq, smat_wide))

#so as cross-loadings of y17 on either factor increase, y3 strength (f1) go up.
#this is after controlling for the fitted loading of the item itself
ggplot(smat_wide, aes(x=y3_gmetric)) + geom_histogram()
ggplot(smat_wide, aes(x=y10_gmetric)) + geom_histogram()
ggplot(smat_wide, aes(y=y10_gmetric, x=y10_fittedloading)) + geom_point() + stat_smooth() #definitely linear, but noisy

summary(lm(y3_gmetric ~ f1c * f2c + f1sq * f2sq + y3_fittedloading, smat_wide))
summary(lm(y10_gmetric ~ f1c * f2c + f1sq * f2sq + y10_fittedloading, smat_wide)) #similar for either comparator

#so: strength for comparators is affected by the target item, but fitted loadings are not
#note that in the simple lm on the aggregate data, the R2 is much higher because we average out *sampling variability* -- we are seeing the average over 100 replications
#thus, if the R2 is higher, it means that the 'signal' of condition predicting metric is stronger by averaging out noise.

summary(lmer(y10_gmetric ~ f1c * f2c + f1sq * f2sq + y10_fittedloading + (1|condition) + (1|graphNum), smat_wide)) #lmer accounting for depending shows minimal differences with lm
summary(lmer(y10_fittedloading ~ f1c * f2c + f1sq * f2sq + (1|condition) + (1|graphNum), smat_wide)) #lmer accounting for depending shows minimal differences with lm: no dependency of loadings on condition
summary(lmer(y3_fittedloading ~ f1c * f2c + f1sq * f2sq + (1|condition) + (1|graphNum), smat_wide)) #lmer accounting for depending shows minimal differences with lm: no dependency of loadings on condition

strength_agg %>% select(y3_gmetric, y10_gmetric, value) %>% cor()

###
##CLOSENESS
#cmat_wide %>% mutate(f1=f1-mean(f1), f2=f2-mean(f2), f1sq=f1^2, f2sq=f2^2) %>% select(f1, f2, f1sq, f2sq) %>% cor() %>% round(2)
#xtabs(~graphNum + poploading, cmat_rand) #check on data structure -- should be 26 per cell. 13 for factor 1 at this population loading, 13 for factor 2 at this population loading (since we have 13^2 cells in the design)
#cmat_wide <- cmat_rand %>% select(condition, graphNum, value, factor, poploading) %>% spread(key=factor, value=poploading) %>% #graphNum contains the identifying column
#cmat_wide <- cmat_rand %>% select(graphNum, value, factor, fittedloading) %>% spread(key=factor, value=fittedloading) %>% #graphNum contains the identifying column
#  mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)
#ggplot(bb_rand, aes(x=f1, y=value)) + geom_point() + stat_smooth()
summary(lm(value ~ f1c + f2c + f1sq + f2sq, cmat_wide)) #so measure becomes hugely predicted by each loading
summary(lmer(value ~ f1c + f2c + f1sq+ f2sq + (1|graphNum), cmat_wide)) #so measure becomes hugely predicted by each loading

#aggregate replications within each combination of f1 and f2loading
closeness_agg <- cmat_wide %>% group_by(f1, f2) %>% summarize_at(vars(value, y10_gmetric, y3_gmetric, y10_fittedloading, y3_fittedloading), funs(mean)) %>% ungroup() %>% mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)
summary(m1 <- lm(value ~ f1c * f2c + f1sq + f2sq, closeness_agg)) #so measure becomes hugely predicted by each loading

#misses the interaction...
apa.reg.table(m1, filename="/Users/mnh5174/TresorSync/Manuscripts/Hallquist Cross-Sectional Networks Critique 2016/sim3 closeness reg.doc")


summary(lm(y10_gmetric ~ f1c * f2c + f1sq + f2sq, closeness_agg)) #so measure becomes hugely predicted by each loading
summary(lm(y10_gmetric ~ f1c * f2c + f1sq + f2sq + y10_fittedloading, closeness_agg)) #so measure becomes hugely predicted by each loading

with(closeness_agg, cor.test(y10_gmetric, value))

ggplot(closeness_agg, aes(x=f1, y=value, color=f2, group=f2)) + geom_line()

#yes, this is interesting....
ggplot(cmat_wide, aes(x=f1, y=value, color=f2, group=f2)) + stat_smooth()

pdf("figures/closeness f1f2 2018.pdf", width=7, height=5)
g3 <- ggplot(closeness_agg, aes(x=f1, y=f2, fill=value)) + geom_tile() +
  #scale_color_gradient() + 
  scale_fill_viridis("Closeness", option="viridis", breaks=c(.06, .08)) +
  xlab("y17 factor 1 loading") + ylab("y17 factor 2 loading") +
  theme_cowplot(font_size=20) + coord_fixed()
#ggtitle("Association of f1 (x) and f2 (y) loadings with strength (color)")
plot(g3)
dev.off()

pdf("figures/closeness f1f2 2018 y10 comparator.pdf", width=7, height=5)
g4 <- ggplot(closeness_agg, aes(x=f1, y=f2, fill=y10_gmetric)) + geom_tile() +
  #scale_color_gradient() + 
  scale_fill_viridis("Closeness", option="viridis", breaks=c(.05, .06)) +
  xlab("y17 factor 1 loading") + ylab("y17 factor 2 loading") +
  theme_cowplot(font_size=20) + coord_fixed()
#ggtitle("Association of f1 (x) and f2 (y) loadings with strength (color)")
plot(g4)
dev.off()

closeness_agg %>% select(y10_gmetric, value) %>% cor()
closeness_agg %>% select(y3_gmetric, value) %>% cor()

#so as cross-loadings of y17 on either factor increase, y3 strength (f1) go up.
#this is after controlling for the fitted loading of the item itself
ggplot(cmat_wide, aes(x=y3_gmetric)) + geom_histogram()
ggplot(cmat_wide, aes(x=y10_gmetric)) + geom_histogram()

summary(lm(y3_gmetric ~ f1c * f2c + f1sq * f2sq + y3_fittedloading, cmat_wide))
summary(lm(y10_gmetric ~ f1c * f2c + f1sq * f2sq + y10_fittedloading, cmat_wide)) #similar for either comparator


##BETWEENNESS
#bmat_wide <- bmat_rand %>% select(graphNum, value, factor, fittedloading) %>% spread(key=factor, value=fittedloading) #use fitted versus population loading
#xtabs(~graphNum + poploading, bmat_rand) #check on data structure -- should be 26 per cell. 13 for factor 1 at this population loading, 13 for factor 2 at this population loading (since we have 13^2 cells in the design)
#need to keep condition and graphNum in the mix to uniquely identify each combination. Otherwise, spread fails when the "value" column is identical
#bmat_wide <- bmat_rand %>% select(condition, graphNum, value, factor, poploading) %>% spread(key=factor, value=poploading) %>% #graphNum contains the identifying column
#  mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)
#ggplot(bmat_wide, aes(x=f1, y=value)) + geom_point() + stat_smooth()
summary(lm(value ~ f1c + f2c + f1sq + f2sq, bmat_wide)) #so measure becomes hugely predicted by each loading
summary(lmer(value ~ f1c + f2c + f1sq + f2sq + (1|condition) + (1|graphNum), bmat_wide)) #so measure becomes hugely predicted by each loading

#aggregate replications within each combination of f1 and f2loading
betweenness_agg <- bmat_wide %>% group_by(f1, f2) %>% summarize_at(vars(value, y10_gmetric, y3_gmetric, y10_fittedloading, y3_fittedloading), funs(mean)) %>% ungroup() %>% mutate(f1c=f1-mean(f1), f2c=f2-mean(f2), f1sq=f1c^2, f2sq=f2c^2)
summary(m1 <- lm(value ~ f1c * f2c, betweenness_agg))
betweenness_agg$vv <- predict(m1)
ggplot(betweenness_agg, aes(x=f1c, y=f2c, fill=vv)) + geom_tile() + coord_fixed()
apa.reg.table(m1, filename="/Users/mnh5174/TresorSync/Manuscripts/Hallquist Cross-Sectional Networks Critique 2016/sim3 betweenness reg.doc")

summary(m2 <- lm(value ~ f1c * f2c + f1sq * f2sq, betweenness_agg)) #so measure becomes hugely predicted by each loading

anova(m1, m2)

summary(m1 <- lm(y10_gmetric ~ f1c * f2c , betweenness_agg)) #+ f1sq * f2sq + y3_fittedloading
apa.reg.table(m1, filename="/Users/mnh5174/TresorSync/Manuscripts/Hallquist Cross-Sectional Networks Critique 2016/sim3 betweenness comparator reg.doc")
summary(m2 <- lm(y10_gmetric ~ f1c * f2c + f1sq * f2sq, betweenness_agg)) #+ f1sq * f2sq + y3_fittedloading
summary(m3 <- lm(y10_gmetric ~ f1c * f2c + f1sq * f2sq + y3_fittedloading, betweenness_agg)) #+ f1sq * f2sq + y3_fittedloading

anova(m1, m2, m3)

summary(lm(y3_gmetric ~ f1c * f2c + f1sq * f2sq + y3_fittedloading, betweenness_agg))
summary(lm(y10_gmetric ~ f1c * f2c + f1sq * f2sq + y10_fittedloading, betweenness_agg)) #similar for either comparator


#yes, this is interesting....
ggplot(bmat_wide, aes(x=f1, y=value, color=f2, group=f2)) + stat_smooth(method="lm")

pdf("figures/betweenness f1f2 2018.pdf", width=7, height=5)
g5 <- ggplot(betweenness_agg, aes(x=f1, y=f2, fill=value)) + geom_tile() +
  #scale_color_gradient() + 
  scale_fill_viridis("Betweenness", option="viridis") +
  xlab("y17 factor 1 loading") + ylab("y17 factor 2 loading") +
  theme_cowplot(font_size=20) + coord_fixed()
#ggtitle("Association of f1 (x) and f2 (y) loadings with strength (color)")
plot(g5)
dev.off()

pdf("figures/betweenness f1f2 2018 y10 comparator.pdf", width=7, height=5)
g6 <- ggplot(betweenness_agg, aes(x=f1, y=f2, fill=y10_gmetric)) + geom_tile() +
  #scale_color_gradient() + 
  scale_fill_viridis("Betweenness", option="viridis") +
  xlab("y17 factor 1 loading") + ylab("y17 factor 2 loading") +
  theme_cowplot(font_size=20) + coord_fixed()
#ggtitle("Association of f1 (x) and f2 (y) loadings with strength (color)")
plot(g6)
dev.off()

#similarity of underlying representation
betweenness_agg %>% select(y10_gmetric, value) %>% cor()
with(betweenness_agg, cor.test(y10_gmetric, value))

#so as cross-loadings of y17 on either factor increase, y3 strength (f1) go up.
#this is after controlling for the fitted loading of the item itself
ggplot(bmat_wide, aes(x=y3_gmetric)) + geom_histogram()
ggplot(bmat_wide, aes(x=y10_gmetric)) + geom_histogram()
summary(lm(y3_gmetric ~ f1c * f2c + f1sq * f2sq + y3_fittedloading, bmat_wide))
summary(lm(y10_gmetric ~ f1c * f2c + f1sq * f2sq + y10_fittedloading, bmat_wide)) #similar for either comparator

##combined figure
pdf("figures/sim3 allmetrics.pdf", width=14, height=17)
plot_grid(g1 + ggtitle("Target"), g2 + ggtitle("Comparator"), 
  g3 + ggtitle("Target"), g4 + ggtitle("Comparator"), 
  g5 + ggtitle("Target"), g6 + ggtitle("Comparator"), ncol=2, align="hv", labels=c("a)", "b)", "c)", "d)", "e)", "f)"), label_size = 22)
dev.off()
