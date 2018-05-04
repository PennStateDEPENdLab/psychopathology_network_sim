library(qgraph)
library(igraph)
library(lavaan)
library(viridis)
library(dplyr)
library(lme4)
library(ggplot2)
library(cowplot)
library(caTools)

setwd(file.path(getMainDir(), "psychopathology_network_sim"))
source("code/psychopathology_sim_functions.R")

load("data/u_vs_g_mlist.RData")

#this is the core setup of the simulation conditions in the mlist
#combine with loading grid to look at association
r2 <- 0.64
f_loading <- sqrt(seq(0, r2, .01))
u_loading <- sqrt(r2 - f_loading^2 + .Machine$double.eps)
f_loading^2 + u_loading^2
loading_grid <- data.frame(condition=1:length(f_loading), f_loading, u_loading, r2_f = f_loading^2, r2_u = u_loading^2)
loading_grid$r2_ratio = with(loading_grid, r2_f/r2_u)
loading_grid$r2_diff = round(with(loading_grid, r2_f - r2_u), 2) #round to avoid tiny diffs
#with(loading_grid[-65,], plot(r2_f, r2_ratio))
#with(loading_grid, plot(r2_f, r2_diff))

#handy bootstrapped mean function
bootf <- function(x, ...) {
  require(rms)
  require(dplyr)
  newstuff <- rbind(smean.cl.boot(x, ...)) %>% data.frame
  return(newstuff)
}



#test equivalence of covariance versus 2-indicator latent
cat(mlist_uniq[[1]]$specification$syntax$fitsyntax)
mm_2indlatent <- cfa(mlist_uniq[[20]]$specification$syntax$fitsyntax, data=mlist_uniq[[65]]$simdata[[30]]) #check the model that was fitted
inspect(mm_2indlatent, 'r2')
#residual correlation syntax
msyn <- '
f1 =~ NA*y1 + NA*y3 + NA*y4 + NA*y5 + NA*y6 + NA*y7 + NA*y8 + NA*y9 + NA*y10 + NA*y2 
f2 =~ NA*y12 + NA*y13 + NA*y14 + NA*y15 + NA*y16 + NA*y17 + NA*y18 + NA*y19 + NA*y20 + NA*y11 
f1 ~~ 1*f1
f1 ~~ 0*f2
f2 ~~ 1*f2
y2 ~~ y11
y1 ~~ start(0.36)*y1
y2 ~~ start(0.36)*y2
y3 ~~ start(0.36)*y3
y4 ~~ start(0.36)*y4
y5 ~~ start(0.36)*y5
y6 ~~ start(0.36)*y6
y7 ~~ start(0.36)*y7
y8 ~~ start(0.36)*y8
y9 ~~ start(0.36)*y9
y10 ~~ start(0.36)*y10
y11 ~~ start(0.36)*y11
y12 ~~ start(0.36)*y12
y13 ~~ start(0.36)*y13
y14 ~~ start(0.36)*y14
y15 ~~ start(0.36)*y15
y16 ~~ start(0.36)*y16
y17 ~~ start(0.36)*y17
y18 ~~ start(0.36)*y18
y19 ~~ start(0.36)*y19
y20 ~~ start(0.36)*y20
'
mm_residcorr <- cfa(msyn, data=mlist_uniq[[65]]$simdata[[5]]) #check the model that was fitted

mean(sapply(mlist_uniq[[1]]$simdata, function(mat) { cor(mat$y2, mat$y11) })) #marginal association of 0.64

mean(sapply(mlist_uniq[[1]]$simdata, function(mat) { summary(lm(y2 ~ y11, mat))$r.squared }))

#I haven't lost my sanity that these are equivalent models
summary(mm_2indlatent,fit.measures=T, standardized=T) #.677^2 = .459 (square factor loading to get correlation)
summary(mm_residcorr,fit.measures=T, standardized=T) #.459 correlation

pdf("figures/sim2_model.pdf", width=20, height=6) #should probably create this in Keynote
semPaths(mm_residcorr, residuals=FALSE,  sizeMan = 4)
dev.off()

#STEP 1: check the relationship between edge strength and factor vs. unique
y2y11_edge <- sapply(mlist_uniq, function(cell) {
    #cell$adjmats$EBICglasso$average["y2", "y11"]
    apply(cell$adjmats$EBICglasso$concat, 1, function(x) { x["y2", "y11"] })
  })

y2y11_ebic_melt <- reshape2::melt(y2y11_edge, varnames=c("replication", "condition"), value.name="edge")
y2y11_ebic_melt$method <- "EBICGLASSO"

y2y11_edge_pcor <- sapply(mlist_uniq, function(cell) {
    #cell$adjmats$pcor$average["y2", "y11"]
    apply(cell$adjmats$pcor$concat, 1, function(x) { x["y2", "y11"] })
  })

y2y11_pcor_melt <- reshape2::melt(y2y11_edge_pcor, varnames=c("replication", "condition"), value.name="edge")
y2y11_pcor_melt$method <- "PCor"

y2y11_all <- rbind(y2y11_ebic_melt, y2y11_pcor_melt)


y2y11_df <- inner_join(y2y11_all, loading_grid, by="condition")

y2y11_byf <- y2y11_df %>% group_by(r2_f, method) %>% dplyr::do(bootf(.$edge, conf.int=0.99)) %>% ungroup() %>% mutate(r2_f = factor(r2_f))
y2y11_byu <- y2y11_df %>% group_by(r2_u, method) %>% dplyr::do(bootf(.$edge, conf.int=0.99)) %>% ungroup() %>% mutate(r2_u = factor(r2_u))

g1 <- ggplot(y2y11_byf, aes(x=r2_f, y=Mean, ymin=Lower, ymax=Upper, color=method, group=method)) + geom_line(size=1) + geom_ribbon(aes(color=NULL), fill="grey60", alpha=0.5) + theme_bw(base_size=14) + xlab("Variance explained by factor") + ylab("Average edge strength between y2 and y11") +
  scale_x_discrete(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) + scale_color_brewer("Edge definition", palette="Dark2") #+ scale_fill_brewer("Edge Method", palette="Dark2")

g2 <- ggplot(y2y11_byu, aes(x=r2_u, y=Mean, ymin=Lower, ymax=Upper, color=method, group=method)) + geom_ribbon(aes(color=NULL), fill="grey88", alpha=1) + geom_line(size=1.3) + theme_bw(base_size=14) + 
  xlab("Variance explained by specific association") + ylab("Average edge strength between y2 and y11") +
  scale_x_discrete(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) + #scale_color_brewer("Edge definition", palette="Dark2") #+ scale_fill_brewer("Edge Method", palette="Dark2")
  scale_color_viridis("Edge definition", discrete=TRUE, begin=0.2, end=0.8, option="cividis")

#pdf("figures/y2y11_association_2panels.pdf", width=10, height=5)
#plot_grid(g1 + theme(legend.position="none", plot.margin=margin(t=30, l=10, b=10, r=15)), g2 + theme(plot.margin=margin(t=30, l=10, b=10, r=15)), align="hv", axis="lt", labels=c("a)", "b)"), rel_widths=c(0.75, 1))
#dev.off()

#on further reflection, because the 2-panel graph is a pure mirror-reverse, just include the residual association plot
pdf("figures/y2y11_association_residonly.pdf", width=7, height=5)
g2 + theme_bw(base_size=16) + theme(axis.title.x = element_text(margin=margin(t=10)), axis.title.y = element_text(margin=margin(r=10)), 
  panel.grid.minor=element_blank(), panel.grid.major=element_line(color="grey94"))
dev.off()


loading_grid$uedge <- colMeans(y2y11_edge)
loading_grid$uedge_pcor <- colMeans(y2y11_edge_pcor)
cor.test(~uedge + r2_f, loading_grid)
cor.test(~uedge + r2_u, loading_grid)

cor.test(~uedge_pcor + r2_f, loading_grid)
cor.test(~uedge_pcor + r2_u, loading_grid)


summary(lm(uedge_pcor ~ r2_f + I(r2_f^2), loading_grid))

summary(m1 <- lm(uedge_pcor ~ r2_u, loading_grid))
summary(m2 <- lm(uedge_pcor ~ r2_u + I(r2_u^2), loading_grid))

anova(m1, m2)

summary(m1 <- lm(uedge ~ r2_u, loading_grid))
summary(m2 <- lm(uedge ~ r2_u + I(r2_u^2), loading_grid))
anova(m1, m2)


plot(~uedge + I(f_loading^2), loading_grid)
plot(~uedge + I(u_loading^2), loading_grid)

plot(~uedge_pcor + I(f_loading^2), loading_grid)
plot(~uedge_pcor + I(u_loading^2), loading_grid)


plot(~uedge + f_loading, loading_grid)
plot(~uedge + u_loading, loading_grid)

####
##step 2: look at graph metrics

nodalstats <- do.call(rbind, lapply(1:length(mlist_uniq), function(i) {
      #updated approach: y2 and y11 are the targeted indicators, but the effects are identical (since the structure is parallel)
      #likewise, all of the other nodes/indicators (y3, y4, etc.) are fungible and their effects look identical, too
      #thus, just focus on y2 and y3 as exemplars for graphs
      #df1 <- filter(mlist_uniq[[i]]$graph_v_factor$EBICglasso$metric_v_loadings, node %in% c("y2", "y3")) %>% mutate(method="EBICGLASSO", condition = i)
      #df2 <- filter(mlist_uniq[[i]]$graph_v_factor$pcor$metric_v_loadings, node %in% c("y2", "y3")) %>% mutate(method="pcor", condition = i)
      df1 <- mlist_uniq[[i]]$graph_v_factor$EBICglasso$metric_v_loadings %>% mutate(method="EBICGLASSO", condition = i)
      df2 <- mlist_uniq[[i]]$graph_v_factor$pcor$metric_v_loadings %>% mutate(method="PCor", condition = i)
      df3 <- mlist_uniq[[i]]$graph_v_factor$cor.shrink$metric_v_loadings %>% mutate(method="Cor.Shrink", condition = i)
      return(rbind(df1, df2, df3))
    }))

#pull in Pearson or cor.shrink correlation for reference? This would require a bigger overhaul since we never defined graphs from marginal association 

#y2 and y11 patterns look identical. Rather than confuse the reader, just show one and mention the other in the caption
nodalstats <- nodalstats %>% inner_join(loading_grid, by="condition")

##ON DECK: look at same graphs for another indicator

#library(forcats)
nodalstats$popr2 <- nodalstats$poploading^2

#tmp <- nodalstats %>% filter(measure=="strength" & method=="EBICGLASSO" & graphNum==1) %>% arrange(condition, graphNum, factor)

#by definition, there is only one unique centrality estimate for the node at every condition that reflects a
#combination of the f (factor) loading and the u (unique) loading. Thus, filter out redundant information.
#there should be one row per condition and metric
nodalstats <- nodalstats %>% group_by(node, measure, method, graphNum, condition) %>% filter(row_number() == 1) %>% ungroup() %>% data.frame() #%>% filter(measure=="strength" & method=="EBICGLASSO") 

#this is now irrelevant since we filter out the redundant info
#nodalstats<- nodalstats %>% mutate(source=recode(factor, "f1"="factor", "f2"="factor", "f3"="residual association"))

#pdf("figures/closeness_vs_loading.pdf", width=12, height=12)
#ggplot(dplyr::filter(nodalstats, measure=="closeness"), aes(y=value, x=popr2)) + geom_point(alpha=0.3) + stat_smooth() +  facet_grid(method~source, scales='free')
#dev.off()

nodalstats_filtered <- nodalstats %>% filter(node %in% c("y2", "y3")) %>% mutate(node=recode(node, y2="Target", y3="Comparator"))

#this is pretty busy with all observations included
pdf("figures/closeness_vs_loading.pdf", width=12, height=12)
g1 <- ggplot(dplyr::filter(nodalstats_filtered, measure=="closeness"), aes(y=value, x=r2_diff)) + geom_point(alpha=0.3) + stat_smooth() + facet_grid(node~method) +
  xlab("Variance explained in y2 (factor - unique)") + ylab("Closeness") + theme_bw(base_size=16)

ggdraw(g1) #+
  #draw_label("More factor", angle = 0, size = 10, x=0.2, y=0.05) +
  #draw_label("More unique", angle = 0, size = 10, x=0.8, y=0.05)
dev.off()

#nrow should go down by a factor of 100 (replications)
nodalstats_repagg <- nodalstats_filtered %>% group_by(node, method, measure, r2_diff) %>% dplyr::do({data.frame(bootf(.$value, conf.int=0.99), r2_diff=.$r2_diff[1], popr2=.$popr2[1])}) %>%
  group_by(node, method, measure) %>% arrange(r2_diff) %>% mutate_at(vars(Lower, Mean, Upper), funs(rm=runmean(., k=5, endrule="mean"))) %>% ungroup()

#create empty plot for re-use
gbase <- ggplot(mapping=aes(y=Mean_rm, ymax=Upper_rm, ymin=Lower_rm, x=r2_diff, color=node)) +
  geom_vline(xintercept=0) +
  geom_line(size=1) + geom_ribbon(alpha=0.2) + facet_wrap(~method, scales="free_y") +
  xlab("Variance explained in Target (factor - unique)") + theme_bw(base_size=18) + 
  scale_color_brewer("Node type", palette="Set1") + theme(
    axis.title.x = element_text(margin=margin(t=12)), axis.title.y = element_text(margin=margin(r=12)),
    legend.key = element_rect(size = 7),
    legend.key.size = unit(1.8, 'lines'),
    panel.spacing=unit(20, "pt"))


g1supp <- gbase %+% dplyr::filter(nodalstats_repagg, measure=="strength") + ylab("Strength")
pdf("figures/strength_vs_loading.pdf", width=16, height=12)
g1 <- gbase %+% dplyr::filter(nodalstats_repagg, measure=="strength" & !method=="Cor.Shrink") + ylab("Strength")
plot(g1)
dev.off()


g2supp <- gbase %+% dplyr::filter(nodalstats_repagg, measure=="closeness") + ylab("Closeness")
pdf("figures/closeness_vs_loading.pdf", width=8, height=6)
g2 <- gbase %+% dplyr::filter(nodalstats_repagg, measure=="closeness" & !method=="Cor.Shrink") + ylab("Closeness")
plot(g2)
dev.off()


g3supp <- gbase %+% dplyr::filter(nodalstats_repagg, measure=="betweenness") + ylab("Betweenness")
pdf("figures/betweenness_vs_loading.pdf", width=12, height=12)
g3 <- gbase %+% dplyr::filter(nodalstats_repagg, measure=="betweenness" & !method=="Cor.Shrink") + ylab("Betweenness") #aes(y=log10(Mean_rm + 1), ymin=log10(Lower_rm + 1), ymax=log10(Upper_rm + 1)) #+ coord_trans(y = "log10")
plot(g3)
dev.off()


legend <- get_legend(g1 + guides(color=guide_legend(title.position="left")) + theme(legend.position="top")) #THESE SEEM TO HAVE NO EFFECT... legend.text=element_text(size=10, margin=margin(l=50, r=50, t=50, b=50)), legend.spacing=unit(100, "pt")
library(cowplot)
#theme(legend.position="none", plot.margin=margin(t=30, l=10, b=10, r=15)), g2 + theme(plot.margin=margin(t=30, l=10, b=10, r=15))

pdf("figures/sim2_nodal_panels.pdf", width=7, height=12)
plot_grid(g1 + xlab("") + theme(legend.position="none"), g2 + xlab("") + theme(legend.position="none"), g3 + theme(legend.position="none"), 
  align="hv", labels=c("a)", "b)", "c)"), axis="ltbr", ncol=1, rel_heights=c(1, 1, 1), label_size=21) #, rel_widths=c(0.75, 1))
#plot_grid(legend, g1 + xlab("") + theme(legend.position="none"), g2 + xlab("") + theme(legend.position="none"), g3 + theme(legend.position="none"), 
#  align="hv", labels=c("", "a)", "b)", "c)"), axis="lt", ncol=1, rel_heights=c(0.3, 1, 1, 1), label_size=21) #, rel_widths=c(0.75, 1))
dev.off()


pdf("figures/sim2_nodal_panels_withcorshrink.pdf", width=9, height=12)
plot_grid(g1supp + xlab("") + theme(legend.position="none"), g2supp + xlab("") + theme(legend.position="none"), g3supp + theme(legend.position="none"), 
  align="hv", labels=c("a)", "b)", "c)"), axis="ltbr", ncol=1, rel_heights=c(1, 1, 1), label_size=21) #, rel_widths=c(0.75, 1))
#plot_grid(legend, g1 + xlab("") + theme(legend.position="none"), g2 + xlab("") + theme(legend.position="none"), g3 + theme(legend.position="none"), 
#  align="hv", labels=c("", "a)", "b)", "c)"), axis="lt", ncol=1, rel_heights=c(0.3, 1, 1, 1), label_size=21) #, rel_widths=c(0.75, 1))
dev.off()

pdf("figures/sim2_nodal_legend.pdf", width=7, height=2)
plot(legend)
dev.off()



nodalstats <- nodalstats %>% mutate(
  r2_diff_sq = r2_diff^2, #don't need to center first because this is already a zero-centered variate (-.64 -- .64)
  r2_diff_cu = r2_diff^3,
  type=case_when(
    node %in% c("y2", "y11") ~ "Target",
    TRUE ~ "Comparator"
  )
)


summary(m0 <- lmer(value ~ 1 + (1|condition:graphNum) + (1|node), filter(nodalstats, measure=="closeness" & method=="EBICGLASSO")))
summary(m1 <- lmer(value ~ r2_diff + type + (1|condition:graphNum), filter(nodalstats, measure=="closeness" & method=="EBICGLASSO")))
summary(m2 <- lmer(value ~ r2_diff + type + (1|condition:graphNum) + (1|node), filter(nodalstats, measure=="closeness" & method=="EBICGLASSO")))
anova(m0, m1, m2) #all variance in node explained by type, it appears (makes sense given design)

summary(m3 <- lmer(value ~ r2_diff + r2_diff_sq + type + (1|condition:graphNum), filter(nodalstats, measure=="closeness" & method=="EBICGLASSO")))
summary(m4 <- lmer(value ~ r2_diff + r2_diff_sq + r2_diff_cu + type + (1|condition:graphNum), filter(nodalstats, measure=="closeness" & method=="EBICGLASSO")))
summary(m4 <- lmer(value ~ r2_diff*type + r2_diff_sq*type + (1|condition:graphNum), filter(nodalstats, measure=="closeness" & method=="EBICGLASSO")))
anova(m1, m3, m4)

library(emmeans)
emtrends(m1, ~type, var="r2_diff")
emtrends(m1, ~type, var="r2_diff_sq")


#STRENGTH
summary(m0 <- lmer(value ~ 1 + (1|condition:graphNum) + (1|node), filter(nodalstats, measure=="strength" & method=="EBICGLASSO")))
summary(m1 <- lmer(value ~ r2_diff + type + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO")))
summary(m2 <- lmer(value ~ r2_diff + type + (1|condition:graphNum) + (1|node), filter(nodalstats, measure=="strength" & method=="EBICGLASSO")))
anova(m0, m1, m2) #all variance in node explained by type, it appears (makes sense given design)

summary(m3 <- lmer(value ~ r2_diff + r2_diff_sq + type + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO")))
summary(m4 <- lmer(value ~ r2_diff + r2_diff_sq + r2_diff_cu + type + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO")))
summary(m5 <- lmer(value ~ r2_diff*type + r2_diff_sq*type + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO")))
summary(m6 <- lmer(value ~ r2_diff*type + r2_diff_sq*type + r2_diff_cu*type + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO")))
anova(m1, m3, m4, m5, m6)

library(emmeans)
emtrends(m6, ~type, var="r2_diff")
emtrends(m6, ~type, var="r2_diff_sq")
emtrends(m6, ~type, var="r2_diff_cu")

#just look at target and comparator separately for a second
summary(m3 <- lmer(value ~ r2_diff + r2_diff_sq + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO" & type=="Target")))
summary(m4 <- lmer(value ~ r2_diff + r2_diff_sq + r2_diff_cu + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO" & type=="Target")))
anova(m3, m4)

summary(m3 <- lmer(value ~ r2_diff + r2_diff_sq + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO" & type=="Comparator")))
summary(m4 <- lmer(value ~ r2_diff + r2_diff_sq + r2_diff_cu + (1|condition:graphNum), filter(nodalstats, measure=="strength" & method=="EBICGLASSO" & type=="Comparator")))
anova(m3, m4)

#too crazy... overpowered. Go to aggregated variant using candidate nodes from graph
nodalstats_repagg <- nodalstats_filtered %>% group_by(node, method, measure, r2_diff) %>% dplyr::do({data.frame(bootf(.$value, conf.int=0.99), r2_diff=.$r2_diff[1], popr2=.$popr2[1])}) %>%
  group_by(node, method, measure) %>% arrange(r2_diff) %>% mutate_at(vars(Lower, Mean, Upper), funs(rm=runmean(., k=5, endrule="mean"))) %>% ungroup() %>%
  mutate(
  r2_diff_sq = r2_diff^2, #don't need to center first because this is already a zero-centered variate (-.64 -- .64)
  r2_diff_cu = r2_diff^3,
  r2_diff_qu = r2_diff^4
)

#STRENGTH GLASSO TARGET
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="strength" & method=="EBICGLASSO" & node=="Target")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="strength" & method=="EBICGLASSO" & node=="Target")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="strength" & method=="EBICGLASSO" & node=="Target")))

anova(m2, m3, m4)

#STRENGTH PCOR TARGET
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="strength" & method=="pcor" & node=="Target")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="strength" & method=="pcor" & node=="Target")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="strength" & method=="pcor" & node=="Target")))

anova(m2, m3, m4)

#STRENGTH GLASSO COMPARATOR
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="strength" & method=="EBICGLASSO" & node=="Comparator")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="strength" & method=="EBICGLASSO" & node=="Comparator")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="strength" & method=="EBICGLASSO" & node=="Comparator")))

anova(m2, m3, m4)

#STRENGTH PCOR COMPARATOR
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="strength" & method=="pcor" & node=="Comparator")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="strength" & method=="pcor" & node=="Comparator")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="strength" & method=="pcor" & node=="Comparator")))

anova(m2, m3, m4)


#closeness GLASSO TARGET
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="closeness" & method=="EBICGLASSO" & node=="Target")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="closeness" & method=="EBICGLASSO" & node=="Target")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="closeness" & method=="EBICGLASSO" & node=="Target")))

anova(m2, m3, m4)

#closeness PCOR TARGET
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="closeness" & method=="pcor" & node=="Target")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="closeness" & method=="pcor" & node=="Target")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="closeness" & method=="pcor" & node=="Target")))

anova(m2, m3, m4)

#closeness GLASSO COMPARATOR
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="closeness" & method=="EBICGLASSO" & node=="Comparator")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="closeness" & method=="EBICGLASSO" & node=="Comparator")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="closeness" & method=="EBICGLASSO" & node=="Comparator")))

anova(m2, m3, m4)

summary(m3 <- lm(Mean ~ r2_diff*node + r2_diff_sq*node , filter(nodalstats_repagg, measure=="closeness" & method=="EBICGLASSO")))
anova(m3)

#closeness PCOR COMPARATOR
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="closeness" & method=="pcor" & node=="Comparator")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="closeness" & method=="pcor" & node=="Comparator")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="closeness" & method=="pcor" & node=="Comparator")))

anova(m2, m3, m4)


#betweenness GLASSO TARGET
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO" & node=="Target")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO" & node=="Target")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO" & node=="Target")))
summary(m5 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu + r2_diff_qu, filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO" & node=="Target")))

anova(m2, m3, m4, m5)

plot(1:65, predict(m5))

#betweenness PCOR TARGET
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="betweenness" & method=="pcor" & node=="Target")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="betweenness" & method=="pcor" & node=="Target")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="betweenness" & method=="pcor" & node=="Target")))
summary(m5 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu + r2_diff_qu, filter(nodalstats_repagg, measure=="betweenness" & method=="pcor" & node=="Target")))

anova(m2, m3, m4, m5)

#betweenness GLASSO COMPARATOR
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO" & node=="Comparator")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO" & node=="Comparator")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO" & node=="Comparator")))

anova(m2, m3, m4)

summary(m3 <- lm(Mean ~ r2_diff*node + r2_diff_sq*node , filter(nodalstats_repagg, measure=="betweenness" & method=="EBICGLASSO")))
anova(m3)

#betweenness PCOR COMPARATOR
summary(m2 <- lm(Mean ~ r2_diff, filter(nodalstats_repagg, measure=="betweenness" & method=="pcor" & node=="Comparator")))
summary(m3 <- lm(Mean ~ r2_diff + r2_diff_sq , filter(nodalstats_repagg, measure=="betweenness" & method=="pcor" & node=="Comparator")))
summary(m4 <- lm(Mean ~ r2_diff + r2_diff_sq + r2_diff_cu, filter(nodalstats_repagg, measure=="betweenness" & method=="pcor" & node=="Comparator")))

anova(m2, m3, m4)




#superseded by above
#
#bwdf <- dplyr::filter(nodalstats, measure=="betweenness") #%>% mutate(logvalue=log10(value + 1.1) + 1)
#
##pdf("figures/betweenness_vs_loading.pdf", width=12, height=12)
##ggplot(bwdf, aes(y=value, x=popr2)) + geom_point(alpha=0.3) + stat_smooth() + facet_wrap(~source, scales='free')
##dev.off()
#
#pdf("figures/betweenness_vs_loading.pdf", width=12, height=12)
#ggplot(bwdf, aes(y=value, x=r2_diff)) + geom_point(alpha=0.3) + stat_smooth() + facet_grid(method~source, scales='free')
#dev.off()
#
#bwdf_agg <- bwdf %>% group_by(method, source, poploading) %>% dplyr::do({data.frame(bootf(.$value, conf.int=0.99), r2_diff=.$r2_diff[1], popr2=.$popr2[1])}) %>% ungroup()
#library(caTools)
#bwdf_agg_rm <- bwdf_agg %>% group_by(method, source) %>% arrange(popr2) %>% mutate_at(vars(Lower, Mean, Upper), funs(runmean(., k=5))) %>% ungroup()
#
#pdf("figures/betweenness_test.pdf", width=12, height=12)
#ggplot(bwdf_agg_rm, aes(y=Mean, x=popr2, ymin=Lower, ymax=Upper)) + geom_line() + geom_ribbon(alpha=0.2) + facet_grid(method~source, scales="free")
#dev.off()
#
##pdf("figures/strength_vs_loading.pdf", width=12, height=12)
##ggplot(dplyr::filter(nodalstats, measure=="strength"), aes(y=value, x=popr2)) + geom_point(alpha=0.3) + stat_smooth() + facet_wrap(~source, scales='free')
##dev.off()
#
#stdf_agg_rm <- dplyr::filter(nodalstats, measure=="strength") %>% group_by(method, source, poploading) %>% 
#  dplyr::do({data.frame(bootf(.$value, conf.int=0.99), .[1,c("r2_diff", "popr2", "condition")])}) %>% ungroup() %>%
#  group_by(method, source) %>% arrange(popr2) %>% mutate_at(vars(Lower, Mean, Upper), funs(runmean(., k=5))) %>% ungroup()  
#
#pdf("figures/strength_vs_loading.pdf", width=12, height=12)
#ggplot(stdf_agg_rm, aes(y=Mean, x=r2_diff, ymin=Lower, ymax=Upper, color=source)) + geom_line() + geom_ribbon(alpha=0.2) + facet_grid(method~., scales='free') +
#  theme_bw(base_size=14)
#dev.off()


summary(m1 <- lmer(value ~ poploading + source + (1|condition:graphNum), filter(nodalstats, measure=="closeness")))
summary(m1 <- lmer(value ~ poploading + source + (1|condition:graphNum), filter(nodalstats, measure=="strength")))
summary(m1 <- lmer(value ~ poploading + source + (1|condition:graphNum), filter(nodalstats, measure=="betweenness")))

nodalstats$r2_diff_sq <- nodalstats$r2_diff^2 
summary(m1 <- lmer(value ~ r2_diff + r2_diff_sq + source + (1|condition:graphNum), filter(nodalstats, measure=="betweenness")))

cm <- lmerCellMeans(m1, yoked=list(list(source="r2_diff", dest="r2_diff_sq", transform=function(x) x^2)))
cm$value <- as.vector(cm$value)
ggplot(cm, aes(x=r2_diff, y=value, color=source, ymin=value-se, ymax=value+se)) + geom_pointrange() + geom_line()



cor.test(~value + poploading, subset(y2eff, factor=="f3" & measure=="betweenness"))
cor.test(~value + poploading, subset(y2eff, factor=="f1" & measure=="betweenness"))

plot(~value + poploading, subset(y2eff, factor=="f3" & measure=="strength"))
plot(~value + poploading, subset(y2eff, factor=="f1" & measure=="strength"))

cor.test(~value + poploading, subset(y2eff, factor=="f3" & measure=="strength"))
cor.test(~value + poploading, subset(y2eff, factor=="f1" & measure=="strength"))

cor.test(~value + poploading, subset(y2eff, factor=="f3" & measure=="closeness"))
cor.test(~value + poploading, subset(y2eff, factor=="f1" & measure=="closeness"))

plot(~value + poploading, subset(y2eff, factor=="f3" & measure=="closeness")) #insanely quadratic
plot(~value + poploading, subset(y2eff, factor=="f1" & measure=="closeness")) #insanely quadratic

qgraph(mlist_uniq[[1]]$adjmats$EBICglasso$average) #0.0 loading on f1, 0.80 loading on u (so y2 and y11 are part of u and not f1)
qgraph(mlist_uniq[[32]]$adjmats$EBICglasso$average) #0.56 loading on f1, 0.57 loading on u (about equal for y2 and y11)
qgraph(mlist_uniq[[65]]$adjmats$EBICglasso$average) #0.8 loading on f1, ~0 loading on u (so y2 and y11 are part of f1 and f2, respectively, and not u) 








#look at how edges for other nodes are affected by the unique y2-y11 edge across loadings
#note that y1 y3 y4 etc. have an equal proportion of variation explained in each case

vv <- do.call(rbind, lapply(1:length(mlist_uniq), function(i) {
      df <- mlist_uniq[[i]]$adjmats$EBICglasso$average
      #df <- df[!rownames(df) %in% c("y2", "y11"), !colnames(df) %in% c("y2", "y11")]
      #df <- df[rownames(df) %in% c("y2", "y11"), colnames(df) %in% c("y2", "y11")]
      #mm <- reshape2::melt(df) %>% mutate(nnum1=as.numeric(substring(as.character(Var1), 2)), nnum2=as.numeric(substring(as.character(Var2), 2))) %>%
      #  filter(nnum1 > nnum2) %>% dplyr::select(-nnum1, -nnum2) %>% mutate(condition=i) #get lower triangle, but keeping node information intact
      
      #just association of targets with comparators
      df <- df[!rownames(df) %in% c("y2", "y11"), colnames(df) %in% c("y2", "y11")] #really need to block this to get on-factor -- otherwise get a mix of on and off
      mm <- reshape2::melt(df) %>% mutate(condition=i) #get lower triangle, but keeping node information intact
      return(mm)
    }))

hist(vv$value)
estrength <- inner_join(vv, loading_grid, by="condition")

summary(m1 <- lmer(value ~ r2_u + (1|condition) , estrength))
summary(m1 <- lm(value ~ r2_u , estrength))

plot(estrength$r2_u, predict(m1))
ggplot(estrength, aes(x=r2_u, y=value)) + geom_point() + stat_smooth()

cm <- lmerCellMeans(m1)
ggplot(cm, aes(x=r2_u, y=value, ymin=value-se, ymax=value+se)) + geom_pointrange()

ggplot(estrength, aes(x=condition, y=estrength, group=condition)) + stat_summary(fun.data = mean_cl_boot, geom = "pointrange")

allf <- mlist_uniq[[65]]$adjmats$EBICglasso$average
allu <- mlist_uniq[[1]]$adjmats$EBICglasso$average
blend <- mlist_uniq[[32]]$adjmats$EBICglasso$average

allf_rmuitems <- allf[!rownames(allf) %in% c("y2", "y11"), !colnames(allf) %in% c("y2", "y11")]
f1_allf <- allf_rmuitems[paste0("y", c(1, 3:10)), paste0("y", c(1, 3:10))]
f1_allf <- f1_allf[lower.tri(f1_allf)]

allu_rmuitems <- allu[!rownames(allu) %in% c("y2", "y11"), !colnames(allu) %in% c("y2", "y11")]
f1_allu <- allu_rmuitems[paste0("y", c(1, 3:10)), paste0("y", c(1, 3:10))]
f1_allu <- f1_allu[lower.tri(f1_allu)]

blend_rmuitems <- blend[!rownames(blend) %in% c("y2", "y11"), !colnames(blend) %in% c("y2", "y11")]
f1_blend <- blend_rmuitems[paste0("y", c(1, 3:10)), paste0("y", c(1, 3:10))]
f1_blend <- f1_blend[lower.tri(f1_blend)]

t.test(f1_blend, f1_allf)
t.test(f1_blend, f1_allu)
t.test(f1_allf, f1_allu)

df <- data.frame(edge=c(f1_blend, f1_allf, f1_allu), 
  condition=rep(c("blend", "allshared", "allunique"), each=length(f1_allf)), subj=1:length(f1_allf))

ggplot(df, aes(x=condition, y=edge)) + geom_boxplot()


#pearson correlation for comparison


library(ez)
aa <- ezANOVA(df, dv=edge, wid=.(subj), within=.(condition))
library(multcomp)
library(lme4)
mm <- lmer(edge ~ condition + (1|subj), df)
summary(glht(mm, linfct = mcp(condition = "Tukey")))
friedman.test(edge ~ condition | subj, df)

allf <- mlist_uniq[[65]]$adjmats$pearson$average
allu <- mlist_uniq[[1]]$adjmats$pearson$average
blend <- mlist_uniq[[32]]$adjmats$pearson$average

allf_rmuitems <- allf[!rownames(allf) %in% c("y2", "y11"), !colnames(allf) %in% c("y2", "y11")]
f1_allf <- allf_rmuitems[paste0("y", c(1, 3:10)), paste0("y", c(1, 3:10))]
f1_allf <- f1_allf[lower.tri(f1_allf)]

allu_rmuitems <- allu[!rownames(allu) %in% c("y2", "y11"), !colnames(allu) %in% c("y2", "y11")]
f1_allu <- allu_rmuitems[paste0("y", c(1, 3:10)), paste0("y", c(1, 3:10))]
f1_allu <- f1_allu[lower.tri(f1_allu)]

blend_rmuitems <- blend[!rownames(blend) %in% c("y2", "y11"), !colnames(blend) %in% c("y2", "y11")]
f1_blend <- blend_rmuitems[paste0("y", c(1, 3:10)), paste0("y", c(1, 3:10))]
f1_blend <- f1_blend[lower.tri(f1_blend)]

t.test(f1_blend, f1_allf)
t.test(f1_blend, f1_allu)
t.test(f1_allf, f1_allu)

df <- data.frame(edge=c(f1_blend, f1_allf, f1_allu), 
  condition=rep(c("blend", "allshared", "allunique"), each=length(f1_allf)), subj=1:length(f1_allf))

ggplot(df, aes(x=condition, y=edge)) + geom_boxplot()

friedman.test(edge ~ condition | subj, df)
ezANOVA(df, dv=edge, wid=.(subj), within=.(condition))