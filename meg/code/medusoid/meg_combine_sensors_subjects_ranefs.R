# Process the RI MEG TF "encoding" model output ("ddf") to get sensor- and subject-wise random coefficients.

plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/"

# first reg_df
regressor <-  "entropy_change"
setwd(paste0(plot_dir, "/", regressor))
ddf <- readRDS(paste0("meg_ddf_wholebrain_", regressor, ".rds"))
# get the necessary estimates, separate ran_vals and ran_coefs into different variables
reg_df <- ddf %>% filter(term == "entropy_change_t" & effect != "fixed") %>%
  select(Time, Freq, group, level, term, effect, estimate) %>%
  pivot_wider(names_from = "effect", values_from = "estimate") 

# get only subject ranefs
subject_df <- reg_df %>% filter(group == "Subject") %>% select(-group) %>% 
  rename(subject_ran_vals = ran_vals, subject_ran_coefs = ran_coefs, subject = level)
# get only sensor ranefs
sensor_df <- reg_df %>% filter(group == "Sensor") %>% select(-group) %>% 
  rename(sensor_ran_vals = ran_vals, sensor_ran_coefs = ran_coefs, sensor = level)


# combine all permutations of sensor and subject, back-calculate fixed effect
fdf <- full_join(subject_df, sensor_df) %>% mutate(combined_ran_vals = subject_ran_vals + sensor_ran_vals,
                                                   fixed_effect = sensor_ran_coefs - sensor_ran_vals,
                                                   combined_effect = fixed_effect + combined_ran_vals)

saveRDS(fdf, "meg_echange_sensor_subject_ranefs.Rds")

regressor <- "reward"
setwd(paste0(plot_dir, "/", regressor))
ddf <- readRDS(paste0("meg_ddf_wholebrain_", regressor, ".rds"))
# get the necessary estimates, separate ran_vals and ran_coefs into different variables
reg_df <- ddf %>% filter(term == "reward_t" & effect != "fixed") %>%
  select(Time, Freq, group, level, term, effect, estimate) %>%
  pivot_wider(names_from = "effect", values_from = "estimate") 

# get only subject ranefs
subject_df <- reg_df %>% filter(group == "Subject") %>% select(-group) %>% 
  rename(subject_ran_vals = ran_vals, subject_ran_coefs = ran_coefs, subject = level)
# get only sensor ranefs
sensor_df <- reg_df %>% filter(group == "Sensor") %>% select(-group) %>% 
  rename(sensor_ran_vals = ran_vals, sensor_ran_coefs = ran_coefs, sensor = level)


# combine all permutations of sensor and subject, back-calculate fixed effect
fdf <- full_join(subject_df, sensor_df) %>% mutate(combined_ran_vals = subject_ran_vals + sensor_ran_vals,
                                                   fixed_effect = sensor_ran_coefs - sensor_ran_vals,
                                                   combined_effect = fixed_effect + combined_ran_vals)
saveRDS(fdf, "meg_reward_sensor_subject_ranefs.Rds")


