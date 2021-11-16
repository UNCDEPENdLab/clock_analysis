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
regressors = c("entropy_change_ri", "abspe_by_rew", "v_max_ri")
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
pos_1 <- all_df %>% filter(term != "abs_pe_sc" & t >= -2 & t <= -1.8 & Freq >= "16.8") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate))

neg_2 <- all_df %>% filter(term != "abs_pe_sc" & t >= .4 & t <= .8 & Freq >= "8.4" & Freq <= "16.8") %>%
  group_by(level) %>%
  summarize(avg=mean(estimate))

cor.test(pos_1$avg, neg_2$avg)
plot(pos_1$avg, neg_2$avg)

################
# same for subjects' betas
# load  models
regressors = c("entropy_change", "v_max", "reward")
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

vm1 <- all_dfs %>% filter(term == "v_max_wi" & t >= 0.4 & t <= 0.8 & Freq >= "8.4" &  Freq <= "16.8") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "vmax_late_beta")

r1 <- all_dfs %>% filter(term == "reward_t" & t >= 0.2 & t <= 0.4 & Freq >= "5" &  Freq <= "8.4") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "reward_early_theta")
r2 <- all_dfs %>% filter(term == "reward_t" & t >= 0.4 & t <= 0.8 &  Freq <= "4.2") %>% # & Freq < "16.8"
  group_by(level) %>%
  summarize(avg=mean(estimate)) %>% mutate(reg_region = "reward_late_delta")
betas <- rbind(ec1, vm1, r1, r2) %>% rename(id = level)
wbetas <- betas %>% pivot_wider(names_from = reg_region, values_from = avg)
setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/")
saveRDS(betas, file = "MEG_betas_echange_vmax_reward_Nov15_2021.RDS")
saveRDS(wbetas, file = "MEG_betas_wide_echange_vmax_reward_Nov15_2021.RDS")

