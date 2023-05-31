##Load Packages and source scripts
pacman::p_load(tidyverse, readr, janitor, ggplot2,wesanderson, cowplot, flextable, plotly, emmeans, lme4, qualtRics, sjPlot, circular, pracma, cowplot)
library(rstan)
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores()-1)
library(brms)
#set home and data directories
#Settings
analysis_dir <- "~/code/clock_analysis/fmri/explore/deepmreye/"
setwd(analysis_dir)

df <- readRDS('deepmreye_explore_n50_processed.rds')
br.formula <- brmsformula(theta_c ~ 1 + (1|id), kappa ~ 1) + von_mises()
ob_R <- brm(br.formula, data = df, chains = 1)
summary(ob_R)

#save(ob, file = paste0(file = paste0(data_dir, 'Data/cannon_brmstidy_4', 'br.ex.Rdata')))
