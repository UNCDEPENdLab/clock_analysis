library(Directional)
library(dplyr)
#library(circmax)
library(circular)
library(brms)
library(ggplot2)
library(tidybayes)
setwd("~/TresorSync/Manuscripts/Paolizzi BPD Uncertainty Masters 2020/Analyses_Feb2022")
home_directory <- "~/TresorSync/Manuscripts/Paolizzi BPD Uncertainty Masters 2020/Analyses_Feb2022"
data_dir <- file.path(home_directory, "Data")
load(file = file.path(data_dir, 'base_effects_data.RData'))
within_effects <- within_effects %>% na.omit()

circ.summary(within_effects$UP, rads = TRUE, fast = FALSE, tol = 1e-07, plot = TRUE)

# requires 30GB RAM?
#xx <- lm.circular(y = within_effects$UP, x = within_effects$PE, type="c-c")

ex_subj <- within_effects %>% filter(ID==100832) %>%
  dplyr::mutate(prior_outcome = dplyr::lag(catch_miss), smallUP = abs(UP) < .2)


#brm_mod2 <- bf(UP ~ PE * cond * prior_outcome + PE*smallUP, kappa ~ 1) + von_mises() #von_mises(link = "tan_half", link_kappa = "identity")
brm_mod2 <- bf(UP ~ PE * cond * prior_outcome, kappa ~ 1) + von_mises() #von_mises(link = "tan_half", link_kappa = "identity")
mm <- brm(brm_mod2, ex_subj)
summary(mm)
plot(mm)

ex_subj %>%
  add_residual_draws(mm) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval()

ex_subj %>%
  add_residual_draws(mm) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()

# this also blows up due to a cartesian product?
# to_plot <- ex_subj %>%
#   add_residual_draws(mm) %>%
#   add_fitted_draws(mm)
# 
# to_plot <- ex_subj %>%
#   add_fitted_draws(mm)

fit_vals <- fitted(mm)[,"Estimate"]
resid_vals <- resid(mm)[,"Estimate"]
plot(fit_vals, resid_vals)

to_check <- ex_subj %>%
  filter(totTrial > 2) %>%
  select(-ID, -fastpredT, -bigPE, -fastpredT_fac, -bigPE_fac, -StaySwitch_fac, -predT_inv, -predT_log) %>%
  mutate(resid = resid_vals, fit = fit_vals, smallUP = abs(UP) < .2)


# this is consistent with the account that the ~0 update trials have some sort of different process or may be invalid in some way
ggplot(to_check, aes(x=fit, y=resid, shape=smallUP, color=predT)) + geom_jitter() + scale_color_viridis_c()

to_check %>% group_by(smallUP) %>% summarise(mean(predT), median(predT))
to_check %>% ggplot(aes(x=predT)) + geom_histogram(bins=12) + facet_wrap(~smallUP)

#try to track down the negative slope line observations
interesting_points <- which(resid_vals > 1 & fit_vals < -0.8)

# all are 0 update, fast RT trials
to_check %>% slice(interesting_points)

#try to track down the negative slope line observations
interesting_points <- which(resid_vals < -0.8 & fit_vals > .8)

to_check %>% slice(interesting_points)

# all are 0 update, fast RT trials
to_check %>% slice(interesting_points)



library(GGally)
pdf("resid_corrs.pdf", width=15, height=15)
ggpairs(to_check)
dev.off()
  
