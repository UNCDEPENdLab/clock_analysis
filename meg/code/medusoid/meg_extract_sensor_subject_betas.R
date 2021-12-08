# correlating sensor ranefs to compare sources, entropy change beta1 vs:
# - late ?abs PE beta
# - early ?abs PE beta
# - Vmax alpha/beta

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
source("~/code/Rhelpers/theme_black.R")

repo_directory <- "~/code/clock_analysis"
data_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/output"
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain"

# load all models
regressors = c("entropy_change_rs","v_max_ri")
all_res <- lapply(regressors, function(rr) {
  cd <- file.path(plot_dir, rr)
  dfile <- file.path(cd, paste0("meg_ddf_wholebrain_", rr, ".rds"))
  ddf <- readRDS(dfile) %>% 
    filter(effect=="ran_vals" & term != "(Intercept)") %>% 
    dplyr::select(t, Freq, level, term, estimate, std.error, conf.low, conf.high)
  return(ddf)
})

all_res[[2]] <- all_res[[2]] %>% filter(term =="abs_pe_sc")
all_res[[1]] <- all_res[[1]] %>% filter(Freq <= 40) %>% droplevels()

all_df <- data.table::rbindlist(all_res)

pe_echange <- all_df %>% filter(term != "v_max_wi" & t >= .4 & t <= .8 & Freq >= "8.4" & Freq < "20") %>%
  select(-std.error, -conf.low, -conf.high) %>%
  pivot_wider(names_from="term", values_from = "estimate") %>%
  group_by(t, Freq) %>%
  do({
    df <- .
    r <- cor(df$entropy_change_t, df$abs_pe_sc)
    data.frame(df[1,], r)
  }) %>% ungroup()

ggplot(pe_echange, aes(x=t, y=Freq, fill=r)) + geom_tile() + scale_fill_viridis_c(begin=0.5)
hist(pe_echange$r)

####
pe_echange <- all_df %>% filter(term != "v_max_wi" & t >= .2 & t <= .4 & Freq >= "4.2" & Freq <= "7") %>%
  select(-std.error, -conf.low, -conf.high) %>%
  pivot_wider(names_from="term", values_from = "estimate") %>%
  group_by(t, Freq) %>%
  do({
    df <- .
    r <- cor(df$entropy_change_t, df$abs_pe_sc)
    data.frame(df[1,], r)
  }) %>% ungroup()

ggplot(pe_echange, aes(x=t, y=Freq, fill=r)) + geom_tile() + scale_fill_viridis_c(begin=0.5)
hist(pe_echange$r)

###
pos_1 <- all_df %>% filter(term == "abs_pe_sc" & t >= .2 & t <= .4 & Freq >= "4.2" & Freq <= "7") %>%
  group_by(level) %>%
  summarize(avg=mean(estimate))

neg_2 <- all_df %>% filter(term == "abs_pe_sc" & t >= .4 & t <= .8 & Freq >= "8.4" & Freq <= "16.8") %>%
  group_by(level) %>%
  summarize(avg=mean(estimate))

cor.test(pos_1$avg, neg_2$avg)
plot(pos_1$avg, neg_2$avg)

################
# same for subjects' betas
# load  models
regressors = c("entropy_change")
# regressors = c("entropy_change", "v_max", "reward")
all_ress <- lapply(regressors, function(rr) {
  cd <- file.path(plot_dir, rr)
  dfile <- file.path(cd, paste0("meg_ddf_wholebrain_", rr, ".rds"))
  ddf <- readRDS(dfile) %>% 
    # filter(effect=="ran_vals" & term != "(Intercept)" & group == "Subject") %>% droplevels() %>% 
    filter(effect=="ran_vals" & group == "Subject") %>% droplevels() %>% 
    dplyr::select(t, Freq, level, term, estimate, std.error, conf.low, conf.high)
  return(ddf)
})

all_dfs <- data.table::rbindlist(all_ress)

# extract entropy change "betas"
# ecdf <- all_dfs %>% filter(term=="entropy_change_t")
ec1 <- all_dfs %>% filter(term == "entropy_change_t" & t >= 0.4 & t <= 0.8 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "entropy_change_late_beta")

ec2 <- all_dfs %>% filter(term == "entropy_change_t" & t >= -0.2 & t <= 0.1 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "entropy_change_early_beta")


# vm1 <- all_dfs %>% filter(term == "v_max_wi" & t >= 0.4 & t <= 0.8 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
#   group_by(level) %>%
#   summarize(avg=mean(estimate)) %>% mutate(reg_region = "vmax_late_beta")
# 
# r1 <- all_dfs %>% filter(term == "reward_t" & t >= 0.2 & t <= 0.4 & Freq >= "5" &  Freq <= "8.4") %>% # & Freq < "16.8"
#   group_by(level) %>%
#   summarize(avg=mean(estimate)) %>% mutate(reg_region = "reward_early_theta")
# r2 <- all_dfs %>% filter(term == "reward_t" & t >= 0.4 & t <= 0.8 &  Freq <= "4.2") %>% # & Freq < "16.8"
#   group_by(level) %>%
#   summarize(avg=mean(estimate)) %>% mutate(reg_region = "reward_late_delta")
betas <- rbind(ec1, ec2) %>% rename(id = level)
# betas <- rbind(ec1, vm1, r1, r2) %>% rename(id = level)
wbetas <- betas %>% pivot_wider(names_from = reg_region, values_from = avg)
setwd("~/code/clock_analysis/meg/data")
# saveRDS(betas, file = "MEG_betas_echange_vmax_reward_Nov15_2021.RDS")
saveRDS(betas, file = "MEG_betas_echange_Nov21_2021.RDS")
saveRDS(wbetas, file = "MEG_betas_wide_echange_Nov21_2021.RDS")


# correlation between fMRI and MEG betas
# get fMRI data
fmri_ec_betas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/fmri_ec_betas.RDS") %>% mutate(id = as.character(id))
multi_betas <- inner_join(wbetas, fmri_ec_betas)

just_betas <- multi_betas %>% select(is.numeric)
cormat <- corr.test(just_betas, method = "pearson")
corrplot(cormat$r, order = "hclust")


# add PE and entropy betas
pebetas <- read.csv("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/pe_max_roi_betas.csv")
pemeta <- read.csv("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/pe_max_cluster_metadata.csv")
pemeta$label <- substr(pemeta$label,22,100)
pemeta_overall <- pemeta[pemeta$l2_contrast == 'overall' & pemeta$l3_contrast == 'Intercept' & pemeta$model == 'Intercept-Age',]
pe <- as_tibble(pebetas[pebetas$l2_contrast == 'overall' & pebetas$l3_contrast == 'Intercept' & pebetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<8 | cluster_number == 10 | cluster_number == 11)
# head(merge(h,meta))
perois_list <- distinct(pemeta_overall[c(1:7, 10:11),c(5,12)])

# inspect distributions
# for some reason, the clusters do not seem to fit the map -- check with MNH (big clusters broken down?) !!
perois <- inner_join(pe,pemeta_overall)
perois$labeled_cluster <- paste(perois$cluster_number,perois$label)
pemeta_overall$labeled_cluster <- paste(pemeta_overall$cluster_number,pemeta_overall$label)
# ggplot(perois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)

pe_labeled <- inner_join(pe,perois_list)
pe_labeled$labeled_cluster <- paste(pe_labeled$cluster_number,pe_labeled$label)
pe_num <- select(pe_labeled,c(1,2,3))
pe_labeled <- select(pe_labeled,c(1,3,5))

pe_wide <- spread(pe_labeled,labeled_cluster,cope_value)
map_df  <- as.tibble(read.csv("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/v_entropy-Intercept_design.txt", sep=""))
map_df$id <- as.character(map_df$ID)
multi_pe_betas <- pe_wide %>% inner_join(map_df %>% select(id, feat_input_id, Age)) %>% inner_join(wbetas)
