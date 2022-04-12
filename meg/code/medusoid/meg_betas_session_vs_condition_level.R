library(tidyverse)
library(lme4)

wbetas_old <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 
  mutate(entropy_change_early_beta_supp_session = -  avg_entropy_change_early_beta,
         entropy_change_late_beta_supp_session = - avg_entropy_change_late_beta,
         omission_early_theta_session = - avg_reward_early_theta) %>% select(
           id, entropy_change_late_beta_supp_session,
                                                                             omission_early_theta_session) 
wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_ec_rewfunc_rt_next_reward_rewfunc_April_5_2022.RDS") %>% 
  mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_ec_rewfunc,
         entropy_change_late_beta_supp = - entropy_change_late_beta_ec_rewfunc,
         rt_shorten_late_beta_supp = - rt_next_late_beta_rt_next,
         omission_early_theta = - reward_early_theta_reward_rewfunc) %>%
  select(id, rewFunc, entropy_change_late_beta_supp, rt_shorten_late_beta_supp, 
         omission_early_theta) 
wwbetas <- wbetas %>% pivot_wider(id_cols = "id", names_from = "rewFunc", 
                                  values_from = c("entropy_change_late_beta_supp", 
                                                  "rt_shorten_late_beta_supp", "omission_early_theta" ))

df <- inner_join(wwbetas, wbetas_old, by = "id")
fbetas <- df %>% select(!id)
cormat <- corr.test(fbetas, method = "pearson")
# note negative correlations between CEVR and other betas
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("session_vs_condition_meg_betas_corr.pdf", height = 10, width = 10)
corrplot(cormat$r, p.mat = cormat$p, order = "hclust", tl.cex = .8, insig = 'blank', method = 'number')
dev.off()



m1 <- lm(entropy_change_late_beta_supp_session ~ entropy_change_late_beta_supp_DEV +
           entropy_change_late_beta_supp_CEV + entropy_change_late_beta_supp_IEV +
           entropy_change_late_beta_supp_CEVR, df)
summary(m1)
Anova(m1)
vif(m1)

m2 <- lm(omission_early_theta_session ~ omission_early_theta_DEV +
           omission_early_theta_CEV + omission_early_theta_IEV +
           omission_early_theta_CEVR, df)
summary(m2)

