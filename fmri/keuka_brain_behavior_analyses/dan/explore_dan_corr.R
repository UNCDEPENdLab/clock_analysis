# explore betwen-region correlations

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(psych)
library(corrplot)
library(corrr)
repo_directory <- "~/code/clock_analysis"

# data & options ----

# data loading options
reprocess = F # otherwise load data from cache
if (!reprocess) {
  tall_only = T # only load wide data (parcels and timepoints as variables)
  wide_only = F
}
# load MEDUSA deconvolved data
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

rt_ts <- clock_comb %>% select(id, run, run_trial, evt_time, label, decon_interp) %>% 
  group_by(id, run, run_trial) %>% arrange(id, run, run_trial) %>%
  pivot_wider(names_from = label, values_from = decon_interp)
str(rt_ts)
labels <- unique(rt_comb$label)
ids <- unique(rt_ts$id)
n = 1
group_cormat <- array(dim = c(length(labels), length(labels), length(ids)))
group_pval <- array(dim = c(length(labels), length(labels), length(ids)))
group_cormat <- list()
for(sub in ids) {
  test <- rt_ts %>% filter(id==sub) %>% ungroup() %>% select(matches(labels))
  cr <- correlate(test, diagonal = 1)
  cr$id <- sub
  group_cormat[[sub]] <- cr
  # group_pval[,,n] <- cormat$p
}
cordf <- do.call(rbind,group_cormat)
allcors <- cordf %>% group_by(term) %>% summarise_at(vars(labels), mean, na.rm = T)
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots')
pdf("clock_corr.pdf", height = 16, width = 16)
allcors %>% network_plot(min_cor = .8,colors = c("white","white","white","white", "white", "red", "blue"))
dev.off()

library(GGally)
ggnet2(allcors)
igraph::graph_from_adjacency_matrix(allcors)

av_pmat <- apply(group_pval, c(1,2), mean)
corrplot(av_cormat, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = av_pmat, sig.level=0.05, insig = "blank")
diag <- nfactors(cormat$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
dan.fa = psych::fa(test, nfactors=2, rotate = "oblimin", fm = "pa")
dan.faba = psych::bassAckward(test, nfactors=2, rotate = "oblimin", fm = "pa")


pe.fa = fa.sort(psych::fa(pejust_rois, nfactors=3))
p