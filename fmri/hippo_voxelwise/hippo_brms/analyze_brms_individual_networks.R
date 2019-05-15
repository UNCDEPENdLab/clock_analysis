#conduct basic analysis of brms-estimated individual networks
library(brms)
library(tidyr)
library(sjstats)
library(igraph)
library(afex)
library(emmeans)

#setwd("/Users/Shared/ics/clock_analysis/fmri/hippo_voxelwise/mplus_var_hippo/")
setwd("/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/hippo_brms")

load("output/m2_brms_rh_4slc.RData")

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

make_graphs <- function(brmsobj, rm_nonsig=TRUE, rm_selfcon=TRUE) {
  require(dplyr)
  require(igraph)
  
  cc <- coef(brmsobj)$id
  cc <- cc[,,!grepl("Intercept", dimnames(cc)[[3]])]
  emp_bayes <- cc[,"Estimate",]
  
  to <- sub("(rh\\d+)_.*", "\\1", colnames(emp_bayes))
  from <- sub("rh\\d+_(rh\\d+).*", "\\1", colnames(emp_bayes))
  
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

gobj <- make_graphs(m2_rh_4slc)

node_df <- gobj$all_nodal

sout_aov <- aov_ez(dv="strength_out", id = "ID", within = "node", data=node_df)
emmeans(sout_aov, ~node)
pairs(emmeans(sout_aov, ~node))

sin_aov <- aov_ez(dv="strength_in", id = "ID", within = "node", data=node_df)
emmeans(sin_aov, ~node)
pairs(emmeans(sin_aov, ~node))

save(gobj, file="m2_brms_rh_4slc_graphmeasures.RData")
