library(brms)
library(tidyr)

library(sjstats)
library(igraph)
#devtools::install_github("garthtarr/edgebundleR")
library(edgebundleR)
library(pacman)
p_load("tidygraph", "ggraph", "scales", "extrafont")

mobj <- "output/rm6_brms_rh_4slc_entropy.RData"
#setwd("/Users/Shared/ics/clock_analysis/fmri/hippo_voxelwise/mplus_var_hippo/")
setwd("/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/hippo_brms")
#load(mobj)

load("output/m2_brms_rh_4slc.RData")

str(rm6_rh_4slc_entropy, max=1)
rm6_rh_4slc_entropy$model
str(rm6_rh_4slc_entropy$ranef)
hdi(rm6_rh_4slc_entropy)
hdi(rm6_rh_4slc_entropy, type="random")
str(rm6_rh_4slc_entropy, max=1)
str(ranef(rm6_rh_4slc_entropy))

vc <- VarCorr(rm6_rh_4slc_entropy)
cmat <- vc$id$cor[,"Estimate",] #just point estimates

#devtools::install_github("garthtarr/edgebundleR")

g <- graph_from_adjacency_matrix(cmat, mode="undirected", weighted = TRUE)

V(g)$name <- sub("v_entropy_wi", "E", V(g)$name, fixed=TRUE)
V(g)$name <- sub("Intercept", "Int", V(g)$name, fixed=TRUE)
gthresh <- delete.edges(g, which(abs(E(g)$weight) < 0.2))
edgebundle(gthresh)

gg <- as_tbl_graph(g)
ggraph(gg, layout='linear', circular=TRUE) + geom_edge_arc(aes(alpha=weight)) +
  geom_node_point(size=3) + 
  theme_graph(foreground="steelblue") +
  coord_fixed()

#just model with no exogenous covariates

vc <- VarCorr(m2_rh_4slc)
cmat <- vc$id$cor[,"Estimate",] #just point estimates

#remove intercepts for the moment
int_names <- (1:ncol(cmat))[grepl("Intercept", dimnames(cmat)[[1]])]
slo_names <- (1:ncol(cmat))[-1*int_names]
smat <- cmat[slo_names, slo_names]

g <- graph_from_adjacency_matrix(smat, mode="undirected", weighted = TRUE, diag = FALSE)

# V(g)$name <- sub("v_entropy_wi", "E", V(g)$name, fixed=TRUE)
# V(g)$name <- sub("Intercept", "Int", V(g)$name, fixed=TRUE)
# gthresh <- delete.edges(g, which(abs(E(g)$weight) < 0.2))
# edgebundle(gthresh)


g <- delete.edges(g, which(abs(E(g)$weight) < 0.2))

gg <- as_tbl_graph(g) %>%
  activate(edges) %>%
  mutate(pos=weight > 0)

ggraph(gg, layout='linear', circular=TRUE) + geom_edge_arc(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue") +
  coord_fixed() + facet_wrap(~pos) +
  scale_edge_color_gradient2("corr", low=muted("blue"), mid="white", high=muted("red"))


ggraph(gg, layout='auto') + geom_edge_fan(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue") +
  coord_fixed() + facet_wrap(~pos) +
  scale_edge_color_gradient2("corr", low=muted("blue"), mid="white", high=muted("red")) +
  ggtitle("Correlations among random slopes for lagged effects")


pos_graph <- ggraph(gg %>% activate(edges) %>% filter(weight > 0), layout='auto') + geom_edge_fan(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue", base_family = "Arial Narrow") +
  coord_fixed() +
  scale_edge_color_gradient2("corr", low=muted("blue"), mid="white", high=muted("red")) +
  ggtitle("Positive")

neg_graph <- ggraph(gg %>% activate(edges) %>% filter(weight < 0), layout='auto') + geom_edge_fan(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue", base_family = "Arial Narrow") +
  coord_fixed() +
  scale_edge_color_gradient2("corr", low=muted("blue"), mid="white", high=muted("red")) +
  ggtitle("Negative")

library(cowplot)
title <- ggdraw() + draw_label("Correlations among lag-1 connectivity random slopes (layout separate)", fontface='bold')

plots <- plot_grid(pos_graph, neg_graph)
pdf("m2_rh_re_corrs.pdf", width=12, height=8)
plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1))
dev.off()



joint_graph <- ggraph(gg, layout='auto') + geom_edge_fan(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue", base_family = "Arial Narrow") +
  coord_fixed() + facet_wrap(~pos) +
  scale_edge_color_gradient2("corr", low=muted("blue"), mid="white", high=muted("red")) +
  ggtitle("Correlations among lag-1 connectivity random slopes (layout together)")

pdf("m2_rh_re_corrs_joint.pdf", width=12, height=8)
plot(joint_graph)
dev.off()

str(coef(m2_rh_4slc))

emp_bayes <- coef(m2_rh_4slc)$id[,"Estimate",]
emp_bayes <- coef(m2_rh_4slc)$id[,"Q97.5",]

sapply(data.frame(emp_bayes), function(x) { 
  (x < 0)/length(x)
})



hypothesis(m2_rh_4slc, "rh6_rh9_l - rh9_rh6_l > 0")

#some random useful plots
stanplot(m2_rh_4slc)

cov2cor(vcov(m2_rh_4slc))
plot(m2_rh_4slc)

#fixed effects
ff <- fixef(m2_rh_4slc)

ff <- ff[!grepl("Intercept", rownames(ff)),]

ests <- ff[,"Estimate"]
to <- sub("(rh\\d+)_.*", "\\1", names(ests))
from <- sub("rh\\d+_(rh\\d+).*", "\\1", names(ests))
sig <- ff[,"Q2.5"] > 0 #no truly sig neg effects in this matrix, so this works

el <- data.frame(from=from, to=to, weight=ests, sig=sig)

g2 <- graph_from_edgelist(as.matrix(el[,1:2]))
E(g2)$weight <- el$weight
E(g2)$sig <- el$sig

gg <- as_tbl_graph(g2) %>%
  activate(edges)

pdf("rh_lag1_connectivity.pdf", width=12, height=12)
ggraph(gg %>% activate(edges) %>% filter(from!=to & sig==TRUE), layout='linear', circular=TRUE) + 
  geom_edge_fan(aes(label = round(weight, 2), color=weight),
                edge_width=1,
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                arrow = arrow(length = unit(3, 'mm')), 
                end_cap = circle(2, 'mm'), check_overlap = TRUE) + 
  
  #geom_edge_arc(aes(color=weight, label=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), nudge_y=.05, nudge_x=-.08, size=8) + 
  theme_graph(foreground="steelblue", base_family = "Arial Narrow") +
  coord_fixed() + 
  scale_edge_color_gradient2("lag1", low=muted("blue"), mid="white", high=muted("red")) +
  ggtitle("Lag-1 directed connectivity")
dev.off()
