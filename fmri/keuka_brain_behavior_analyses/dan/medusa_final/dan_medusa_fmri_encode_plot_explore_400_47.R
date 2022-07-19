# plot Explore MEDUSA replication results, all models, for 47 DAN nodes in Schaefer 400-17 parcellation

library(tidyverse)
library(lme4)
library(data.table)
library(readxl)
library(fmri.pipeline)

if (Sys.getenv("USER")=="alexdombrovski" | Sys.getenv("USER")=="Alex") {
  repo_directory <- "~/code/clock_analysis"} else {
    repo_directory <- "~/Data_Analysis/clock_analysis"}
out_dir <- file.path(paste0(repo_directory, "/fmri/keuka_brain_behavior_analyses/dan/medusa_final/"))

# read in Explore MEDUSA results
setwd(out_dir)
explore_ddf <- readRDS("rt_parcel_group_all_models_explore_400_47_encode_July19_2022.rds")

# FDR correction
explore_ddf <- explore_ddf %>% dplyr::filter(evt_time <= 5) %>% 
  filter(effect=="fixed") %>%
  group_by(term, model_name) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

# plot all models
plot_medusa(explore_ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="parcel_group", facet_by="side",
            out_dir=file.path(out_dir, "explore_400_47"),  p.value="padj_BY_term")
