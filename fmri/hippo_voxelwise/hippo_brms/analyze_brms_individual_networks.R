#conduct basic analysis of brms-estimated individual networks
library(brms)
library(tidyr)
library(sjstats)
library(igraph)
library(afex)
library(emmeans)
library(ggplot2)
library(cowplot)
library(tidygraph)
library(ggraph)
library(scales)


#setwd("/Users/Shared/ics/clock_analysis/fmri/hippo_voxelwise/mplus_var_hippo/")
#setwd("/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/hippo_brms")
setwd("~/Data_Analysis/clock_analysis/fmri/hippo_voxelwise")

#load("output/m2_brms_rh_4slc.RData")
load("output/lm1_brms_lh_6slc.RData")
load("output/rm1_brms_rh_6slc.RData")

#construct empirical bayes networks
#nodal graph measure worker function
gmeasures <- function(g) {
  gpos <- g
  if (any(E(g)$weight < 0)) {
    message("For ID: ", g$ID, ", removing negative weights in graph: ", paste(E(g)[which(E(g)$weight < 0)], collapse=", "))
    gpos <- delete.edges(g, which(E(g)$weight < 0)) 
  }
  
  evc <- evcent(g)$vector
  sin <- strength(g, mode = "in")
  sout <- strength(g, mode = "out")
  #cc <- closeness(gpos) #too finicky due to disconnected subgraphs
  bet <- betweenness(gpos, normalize=TRUE)*100
  locclust <- transitivity(g, type="barrat")
  pr <- page.rank(g, algo="prpack")$vector
  d <- graph.density(g) #add observed density (graph-level) to nodal stats for further analysis
  if (!is.null(E(g)$weight)) {
    meanW <- mean(E(g)$weight, na.rm=TRUE)
  } else {
    meanW <- NA
  }
  
  #return a 264 x nmetrics matrix for analysis/reduction            
  data.frame(ID=g$ID, node=V(g)$name, density=d, evcent=evc, betweenness=bet, #closeness=cc, 
             locclust=locclust, pagerank=pr, strength_out=sout, strength_in=sin, meanW=meanW)
}

make_graphs <- function(brmsobj, roi_name="rh", rm_nonsig=TRUE, rm_selfcon=TRUE) {
  require(dplyr)
  require(igraph)
  
  cc <- coef(brmsobj)$id
  cc <- cc[,,!grepl("Intercept", dimnames(cc)[[3]])]
  emp_bayes <- cc[,"Estimate",]
  
  to <- sub(paste0("(", roi_name, "\\d+)_.*"), "\\1", colnames(emp_bayes))
  from <- sub(paste0(roi_name, "\\d+_(", roi_name, "\\d+).*"), "\\1", colnames(emp_bayes))
  
  glist <- list()
  for (i in 1:nrow(emp_bayes)) {
    sig <- sapply(1:ncol(emp_bayes), function(col) {
      if (emp_bayes[i,col] < 0 && cc[i,"Q97.5",col] < 0) {
        sig <- TRUE #sig negative
      } else if (emp_bayes[i,col] > 0 && cc[i,"Q2.5",col] > 0) {
        sig <- TRUE #sig positive
      } else {
        sig <- FALSE
      }
    })
    
    el <- data.frame(from=from, to=to, weight=emp_bayes[i,], sig=sig)
    if (rm_selfcon) { el <- el %>% filter(from != to) }
    gg <- graph_from_edgelist(as.matrix(el[,1:2]))
    E(gg)$weight <- el$weight
    E(gg)$sig <- el$sig
    if (rm_nonsig) { gg <- delete.edges(gg, which(E(gg)$sig == FALSE)) }
    gg$ID <- rownames(emp_bayes)[i]
    glist[[i]] <- gg
  }
  
  all_nodal <- do.call(rbind, lapply(glist, gmeasures))
  return(list(glist=glist, all_nodal=all_nodal))
  
  
}

#gobj <- make_graphs(m2_rh_4slc)

gobj_lh <- make_graphs(lm1_lh_6slc, roi_name="lh")
gobj_rh <- make_graphs(rm1_rh_6slc, roi_name="rh")

node_df_lh <- gobj_lh$all_nodal
node_df_rh <- gobj_rh$all_nodal

#basic strength comparisons
#left
sout_aov <- aov_ez(dv="strength_out", id = "ID", within = "node", data=node_df_lh)
emmeans(sout_aov, ~node)
pairs(emmeans(sout_aov, ~node))

sin_aov <- aov_ez(dv="strength_in", id = "ID", within = "node", data=node_df_lh)
emmeans(sin_aov, ~node)
pairs(emmeans(sin_aov, ~node))

#right
sout_aov <- aov_ez(dv="strength_out", id = "ID", within = "node", data=node_df_rh)
emmeans(sout_aov, ~node)
pairs(emmeans(sout_aov, ~node))

sin_aov <- aov_ez(dv="strength_in", id = "ID", within = "node", data=node_df_rh)
emmeans(sin_aov, ~node)
pairs(emmeans(sin_aov, ~node))


#plot edge weights
#keep non-sig edges for edge regressions
gobj_all_lh <- make_graphs(lm1_lh_6slc, roi_name="lh", rm_nonsig=FALSE)
gobj_all_rh <- make_graphs(rm1_rh_6slc, roi_name="rh", rm_nonsig=FALSE)

#plot edge strength by hemisphere
####LEFT
edge_df_lh <- do.call(rbind, lapply(gobj_all_lh$glist, function(g) {
  df <- get.data.frame(g)
  df$ID <- g$ID
  return(df)
} )) %>%
  mutate(from=sub("([A-z]{2})([0-9]{1})$", "\\10\\2", from), to=sub("([A-z]{2})([0-9]{1})$", "\\10\\2", to)) #zero pad single-digit cases for a->p order in plot

edist <- ggplot(edge_df_lh, aes(x=weight)) + geom_density() + facet_grid(rows=vars(from), cols=vars(to)) + geom_vline(xintercept=0, color='orange') +
  theme_bw(base_size=12) + theme(panel.spacing = unit(1.5, "lines"), plot.margin=margin(1.5, 1.5, 1.5, 1.5, 'cm'))

pdf("lh_edistributions.pdf", width=12, height=12)
toplot <- ggdraw(edist) + draw_label("From", angle=-90, size=18, x=0.98, y=0.5) +
  draw_label("To", angle=0, size=18, x=0.5, y=0.98)
plot(toplot)
dev.off()

edge_df_lh$eid <- paste0(edge_df_lh$from, edge_df_lh$to)
e_test <- lmer(weight ~ eid + (1|ID), edge_df_lh)
summary(e_test)
emmeans(e_test, ~eid)

#break out from -> to strength differences by slice
e_test <- lmer(weight ~ from*to + (1|ID), edge_df_lh, REML=FALSE)
summary(e_test)
emmeans(e_test, ~to | from)
pairs(emmeans(e_test, ~to | from))

df <- summary(emmeans(e_test, ~to * from))
gg_lh <- graph_from_edgelist(as.matrix(df[,c("from", "to")]))
E(gg_lh)$weight <- df$emmean

lmer_lh_plot <- as_tbl_graph(gg_lh) %>%
  activate(edges) %>%
  filter(!is.na(weight))

ggraph(lmer_lh_plot, layout='linear', circular=TRUE) + geom_edge_arc(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue") +
  coord_fixed() + 
  scale_edge_color_gradient2("weight", low=muted("blue"), mid="white", high=muted("red"))


png("lh_brms_group_weights.png", width=12, height=12, units="in", res=300)
ggraph(lmer_lh_plot, layout='auto') + 
  geom_edge_fan(edge_width=1.5, aes(color=weight, start_cap = label_rect(node1.name), end_cap = label_rect(node2.name)), arrow = arrow(type = "open", length = unit(3, 'mm'))) +
  #geom_edge_fan(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue") +
  #coord_fixed() +
  scale_edge_color_gradient2("weight", low=muted("blue"), mid="white", high=muted("red")) +
  ggtitle("LH average edge weight from BRMS")
dev.off()


###RIGHT
edge_df_rh <- do.call(rbind, lapply(gobj_all_rh$glist, function(g) {
  df <- get.data.frame(g)
  df$ID <- g$ID
  return(df)
} )) %>%
  mutate(from=sub("([A-z]{2})([0-9]{1})$", "\\10\\2", from), to=sub("([A-z]{2})([0-9]{1})$", "\\10\\2", to)) #zero pad single-digit cases for a->p order in plot

edist <- ggplot(edge_df_rh, aes(x=weight)) + geom_density() + facet_grid(rows=vars(from), cols=vars(to)) + geom_vline(xintercept=0, color='orange') +
  theme_bw(base_size=12) + theme(panel.spacing = unit(1.5, "lines"), plot.margin=margin(1.5, 1.5, 1.5, 1.5, 'cm'))

pdf("rh_edistributions.pdf", width=12, height=12)
toplot <- ggdraw(edist) + draw_label("From", angle=-90, size=18, x=0.98, y=0.5) +
  draw_label("To", angle=0, size=18, x=0.5, y=0.98)
plot(toplot)
dev.off()

edge_df_rh$eid <- paste0(edge_df_rh$from, edge_df_rh$to)
e_test <- lmer(weight ~ eid + (1|ID), edge_df_rh)
summary(e_test)
emmeans(e_test, ~eid)

#break out from -> to strength differences by slice
e_test <- lmer(weight ~ from*to + (1|ID), edge_df_rh, REML=FALSE)
summary(e_test)
emmeans(e_test, ~to | from)
pairs(emmeans(e_test, ~to | from))

df <- summary(emmeans(e_test, ~to * from))
gg_rh <- graph_from_edgelist(as.matrix(df[,c("from", "to")]))
E(gg_rh)$weight <- df$emmean

lmer_rh_plot <- as_tbl_graph(gg_rh) %>%
  activate(edges) %>%
  filter(!is.na(weight))

ggraph(lmer_rh_plot, layout='linear', circular=TRUE) + geom_edge_arc(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue") +
  coord_fixed() + 
  scale_edge_color_gradient2("weight", low=muted("blue"), mid="white", high=muted("red"))


png("rh_brms_group_weights.png", width=12, height=12, units="in", res=300)
ggraph(lmer_rh_plot, layout='auto') + 
  geom_edge_fan(edge_width=1.5, aes(color=weight, start_cap = label_rect(node1.name), end_cap = label_rect(node2.name)), arrow = arrow(type = "open", length = unit(3, 'mm'))) +
  #geom_edge_fan(aes(color=weight)) +
  geom_node_point(aes(label=name), size=3) + 
  geom_node_text(aes(label=name), repel=TRUE) + 
  theme_graph(foreground="steelblue") +
  #coord_fixed() +
  scale_edge_color_gradient2("weight", low=muted("blue"), mid="white", high=muted("red")) +
  ggtitle("RH average edge weight from BRMS")
dev.off()


rh_graph <- as_tbl_graph(gg_rh) %>% activate(edges) %>% mutate(side="Right") %>%
  activate(nodes) %>% mutate(name=sub("^.", "", name))
lh_graph <- as_tbl_graph(gg_lh) %>% activate(edges) %>% mutate(side="Left") %>%
  activate(nodes) %>% mutate(name=sub("^.", "", name))
comb_graph <- graph_join(rh_graph, lh_graph) %>%
  activate(edges) %>%
  filter(!is.na(weight) & abs(weight) > 0.015)

pdf("Brms_6slc_combined_group_graphs_filtered.pdf", width=15, height=8)
ggraph(comb_graph, layout='auto') + #geom_edge_arc(aes(color=weight)) +
  geom_edge_arc(edge_width=1.5, aes(color=weight, start_cap = label_rect(node1.name), end_cap = label_rect(node2.name)), 
                arrow = arrow(type = "open", length = unit(3, 'mm'))) +
  geom_node_point(size=3) + 
  geom_node_text(aes(label=name), color="black", size=9, repel=TRUE) + 
  theme_graph(foreground="steelblue") +
  facet_edges(~side) +
  scale_edge_color_gradient2("weight", low=muted("blue"), mid="white", high=muted("red")) +
  #scale_edge_color_gradient2("weight", low="#87CEFA", mid="gray60", high="#EEE8AA") +
  theme_graph(background="white", caption_colour="gray60", strip_text_size = 20,
              strip_text_colour="gray30", base_size=15)
dev.off()

#e_test_brms <- brm(weight ~ from*to + (1|ID), edge_df_rh) # brms of brms? only in your dreams
#summary(e_test_brms)

#save(gobj_lh, file="output/lm1_lh_6slc_brms_graphmeasures.RData")
#save(gobj_rh, file="output/rm1_rh_6slc_brms_graphmeasures.RData")

save(gobj_lh, gobj_rh, gobj_all_lh, gobj_all_rh, edge_df_lh, edge_df_rh, file="output/brms_6slc_graphs.RData")

edge_df_rh_long <- edge_df_rh %>% dplyr::rename(outgoing=from, incoming=to) %>% gather(key="direction", value="slc", outgoing, incoming)
pdf("brms_6slc_rh.pdf", width=8, height=6)
ggplot(edge_df_rh_long, aes(x=slc, y=weight, color=direction)) + geom_point(position=position_dodge(width=0.4), alpha=0.05) + stat_summary(fun.data = "mean_cl_boot", geom="crossbar") +
  ylim(c(-0.1, 0.18))
dev.off()


edge_df_lh_long <- edge_df_lh %>% dplyr::rename(outgoing=from, incoming=to) %>% gather(key="direction", value="slc", outgoing, incoming)
pdf("brms_6slc_lh.pdf", width=8, height=6)
ggplot(edge_df_lh_long, aes(x=slc, y=weight, color=direction)) + geom_point(position=position_dodge(width=0.4), alpha=0.05) + stat_summary(fun.data = "mean_cl_boot", geom="crossbar") +
  ylim(c(-0.1, 0.18))
dev.off()
