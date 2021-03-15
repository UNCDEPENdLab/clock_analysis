# plot decon time courses by compresson level

library(modelr)
library(tidyverse)
library(lme4)
library(afex)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)
library(viridis)
library(psych)
library(corrplot)
repo_directory <- "~/code/clock_analysis"

# data & options ----

# data loading options
reprocess = F # otherwise load data from cache
if (!reprocess) {
  wide_only = F # only load wide data (parcels and timepoints as variables)
  tall_only = T
}
replicate_compression = F
if(replicate_compression) {reprocess = T}
# load MEDUSA deconvolved data
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

# add run-wise compression
comp <- as_tibble(read_csv("~/Box/SCEPTIC_fMRI/dan_medusa/run_level_compression.csv.gz"))

# merge
clock_comb <- merge(clock_comb, comp)
rt_comb <- merge(rt_comb, comp)
describe(comp$acom_mean)


# within-trial plots
setwd("~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots")
pdf("dan_deconBYcompressonBYregion.pdf", height = 16, width = 20)
ggplot(clock_comb %>% filter(!is.na(v_entropy_wi)), aes(evt_time, decon_interp, color = acom_mean >.66, lty = v_entropy_wi > 0)) + geom_smooth(method = "loess", se = F) + facet_wrap(~label)
dev.off()

