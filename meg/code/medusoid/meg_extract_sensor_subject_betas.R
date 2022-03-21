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
library(psych)
library(corrplot)
source("~/code/Rhelpers/theme_black.R")

repo_directory <- "~/code/clock_analysis"
data_dir <- "/Users/alexdombrovski/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/plots/wholebrain/output"
plot_dir <- "/Users/alexdombrovski/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/plots/wholebrain"
stat_summaries <- T # get mean statistics for every signal

sensors = F # extract sensor random slopes



fmri = T # whether to merge and correlate with fMRI betas
# select_ec_sensors = T
################
# subjects' betas
# load  models
regressors = c("entropy_change", "v_max", "reward", "signed_pe_rs")
# regressors = c("entropy_change", "entropy_change_ec_sensors", "v_max", "reward", "abspe_ec_sensors", "abs_pe")
# regressors = c("entropy_change", "v_max", "reward")
all_regs <- lapply(regressors, function(rr) {
  cd <- file.path(plot_dir, rr)
  if (str_detect(rr, "ec_sensors")) {
    dfile <- file.path(cd, paste0("meg_mixed_by_tf_ddf_ec_sensors_", rr, ".rds"))
  }
  dfile <- file.path(cd, paste0("meg_ddf_wholebrain_", rr, ".rds"))
  ddf <- readRDS(dfile) %>% 
    # filter(effect=="ran_vals" & term != "(Intercept)" & group == "Subject") %>% droplevels() %>% 
    filter(effect=="ran_vals" & group == "Subject") %>% droplevels() %>% 
    dplyr::select(t, Freq, level, term, estimate, std.error, conf.low, conf.high) %>% mutate(regressor = rr)
  return(ddf)
})

all_dfs <- data.table::rbindlist(all_regs)


# extract entropy change "betas"
# ecdf <- all_dfs %>% filter(term=="entropy_change_t")
# late beta
ec1 <- all_dfs %>% filter(term == "entropy_change_t" & t >= 0.5 & t <= 0.8 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
  group_by(level, regressor) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "entropy_change_late_beta")
# early beta (rebound)
ec2 <- all_dfs %>% filter(term == "entropy_change_t" & t >= -0.2 & t <= 0.1 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
  group_by(level, regressor) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "entropy_change_early_beta")

# go up to 20 Hz, arbitrary boundary of beta1
ec1_20 <- all_dfs %>% filter(term == "entropy_change_t" & t >= 0.5 & t <= 0.8 & Freq >= "8.4" &  Freq <= "20") %>% # & Freq < "16.8"
  group_by(level, regressor) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "entropy_change_late_beta_20hz")

ec2_20 <- all_dfs %>% filter(term == "entropy_change_t" & t >= -0.2 & t <= 0.1 & Freq >= "8.4" &  Freq <= "20") %>% # & Freq < "16.8"
  group_by(level, regressor) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "entropy_change_early_beta_20hz")


ecs <- rbind(ec1, ec2, ec1_20, ec2_20) %>% rename(id = level) %>% ungroup() %>% pivot_wider(names_from = c(reg_region, regressor), values_from = avg)


vm1 <- all_dfs %>% filter(term == "v_max_wi" & t >= 0.4 & t <= 0.8 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "vmax_late_beta")

r1 <- all_dfs %>% filter(term == "reward_t" & t >= 0.2 & t <= 0.4 & Freq >= "5" &  Freq <= "8.4") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "reward_early_theta")
r2 <- all_dfs %>% filter(term == "reward_t" & t >= 0.4 & t <= 0.8 &  Freq <= "4.2") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "reward_late_delta")

pe1 <- all_dfs %>% filter(term == "pe_max" & regressor == "signed_pe_rs" & t >= 0.2 & t <= 0.4 & Freq >= "5" &  Freq <= "8.4") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "pe_early_theta")


pe2 <- all_dfs %>% filter(term == "pe_max" & regressor == "signed_pe_rs" & t >= 0.5 & t <= 0.8 &  Freq >= "8.4" &  Freq <= "20") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "pe_late_beta")


# a1 <- all_dfs %>% filter(term == "scale(abs_pe)" & regressor == "abs_pe" & t >= 0.2 & t <= 0.4 & Freq >= "5" &  Freq <= "8.4") %>% # & Freq < "16.8"
#   group_by(level) %>%
#   summarize(avg=mean(estimate)) %>% mutate(reg_region = "abspe_early_theta")
# 
# 
# a2 <- all_dfs %>% filter(term == "scale(abs_pe)" & regressor == "abs_pe" & t >= 0.4 & t <= 0.75 &  Freq >= "8.4" &  Freq <= "20") %>% # & Freq < "16.8"
#   group_by(level) %>%
#   summarize(avg=mean(estimate)) %>% mutate(reg_region = "abspe_late_beta")

# betas <- rbind(ec1, ec2) %>% rename(id = level) %>% ungroup()
wholebrain_betas <- rbind(vm1, r1, r2, pe1, pe2) %>% rename(id = level) %>% pivot_wider(names_from = c(reg_region), values_from = avg) %>% ungroup()
wbetas <- merge(ecs, wholebrain_betas)
just_meg_betas <- wbetas %>% select(where(is.numeric))
cormat <- corr.test(just_meg_betas, method = "pearson")
# inspect corr structure to make sure it is not garbage
corrplot(cormat$r, p.mat = cormat$p, order = "hclust", tl.cex = .8, insig = 'blank', method = 'number')
setwd("~/code/clock_analysis/meg/data")
# saveRDS(betas, file = "MEG_betas_echange_vmax_reward_Nov15_2021.RDS")


# saveRDS(betas, file = paste("MEG_betas", paste(regressors, collapse="_"), "Dec15_2021.RDS", sep = "_"))
saveRDS(wbetas, file = paste("MEG_betas", paste(regressors, collapse="_"), "Mar_14_2022.RDS", sep = "_"))

if (stat_summaries) {

all_stats <- lapply(regressors, function(rr) {
  cd <- file.path(plot_dir, rr)
  dfile <- file.path(cd, paste0("meg_ddf_wholebrain_", rr, ".rds"))
  sdf <- readRDS(dfile) %>% 
    # filter(effect=="ran_vals" & term != "(Intercept)" & group == "Subject") %>% droplevels() %>% 
    filter(effect=="fixed") %>% droplevels() %>% 
    dplyr::select(t, Freq, term, statistic, estimate, std.error, conf.low, conf.high, rhs) %>% mutate(regressor = rr)
  return(sdf)
})

stat_df <- data.table::rbindlist(all_stats) %>% filter(
  (regressor == "entropy_change" & term == "entropy_change_t") |
    (regressor == "v_max" & term == "v_max_wi") | 
    (regressor == "reward" & term == "reward_t") | 
    (regressor == "signed_pe_rs" & term == "pe_max")  
)
}
setwd(plot_dir)
pdf("Random_slope_TF_plots_echange_wi_pe_max_reward_v_max_wi_t_3.pdf", height = 4, width = 6)
ggplot(stat_df %>% filter(abs(statistic) > 3), aes(t, Freq, fill = statistic)) + geom_tile() + facet_wrap(~term) + 
  scale_fill_viridis_c("Statistic,\nthresholded\nat 3") + geom_vline(xintercept = 0, linetype="dotted", 
                                                                   color = "red", size=1.5)
dev.off()

# get mean statistics

sec1 <- stat_df %>% filter(term == "entropy_change_t" & t >= 0.5 & t <= 0.8 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
    summarize(stat=mean(statistic)) %>% mutate(reg_region = "entropy_change_late_beta")
# early beta (rebound)
sec2 <- stat_df %>% filter(term == "entropy_change_t" & t >= -0.2 & t <= 0.1 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
    summarize(stat=mean(statistic)) %>% mutate(reg_region = "entropy_change_early_beta")

# go up to 20 Hz, arbitrary boundary of beta1
sec1_20 <- stat_df %>% filter(term == "entropy_change_t" & t >= 0.5 & t <= 0.8 & Freq >= "8.4" &  Freq <= "20") %>% # & Freq < "16.8"
  summarize(stat=mean(statistic)) %>% mutate(reg_region = "entropy_change_late_beta_20hz")

sec2_20 <- stat_df %>% filter(term == "entropy_change_t" & t >= -0.2 & t <= 0.1 & Freq >= "8.4" &  Freq <= "20") %>% # & Freq < "16.8"
  summarize(stat=mean(statistic)) %>% mutate(reg_region = "entropy_change_early_beta_20hz")

svm1 <- stat_df %>% filter(term == "v_max_wi" & t >= 0.4 & t <= 0.8 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
  summarize(stat=mean(statistic)) %>% mutate(reg_region = "vmax_late_beta")

sr1 <- stat_df %>% filter(term == "reward_t" & t >= 0.2 & t <= 0.4 & Freq >= "5" &  Freq <= "8.4") %>% # & Freq < "16.8"
  summarize(stat=mean(statistic)) %>% mutate(reg_region = "reward_early_theta")
sr2 <- stat_df %>% filter(term == "reward_t" & t >= 0.4 & t <= 0.8 &  Freq <= "4.2") %>% # & Freq < "16.8"
  summarize(stat=mean(statistic)) %>% mutate(reg_region = "reward_late_delta")
spe1 <- stat_df %>% filter(term == "pe_max" & regressor == "signed_pe_rs" & t >= 0.2 & t <= 0.4 & Freq >= "5" &  Freq <= "8.4") %>% # & Freq < "16.8"
  summarize(stat=mean(statistic)) %>% mutate(reg_region = "pe_early_theta")
spe2 <- stat_df %>% filter(term == "pe_max" & regressor == "signed_pe_rs" & t >= 0.5 & t <= 0.8 &  Freq >= "8.4" &  Freq <= "20") %>% # & Freq < "16.8"
  summarize(stat=mean(statistic)) %>% mutate(reg_region = "pe_late_beta")

signal_stats <- rbind(sec1, sec2, sec1_20, sec2_20, svm1, sr1, sr2, spe1, spe2) 

signal_stats_wide <- signal_stats %>% pivot_wider(names_from = c(reg_region), values_from = stat)

pdf("meg_signals_stat_comparison.pdf", height = 4, width = 4)
ggplot(signal_stats, aes(reg_region, abs(stat))) + geom_point() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
# correlation between fMRI and MEG betas
# get fMRI data

if (fmri) {
fmri_ec_betas <- readRDS("~/code/clock_analysis/meg/data/fmri_ec_betas.RDS") %>% mutate(id = as.character(id))
multi_betas <- inner_join(wbetas, fmri_ec_betas)

just_betas <- multi_betas %>% select(is.numeric)
cormat <- corr.test(just_betas, method = "pearson")
corrplot(cormat$r, p.mat = cormat$p, order = "hclust", tl.cex = .5, method = 'number', sig.level = .1, insig = 'blank')
}
if (sensors) {
  # load all models
  regressors = c("entropy_change","v_max_ri")
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
}
