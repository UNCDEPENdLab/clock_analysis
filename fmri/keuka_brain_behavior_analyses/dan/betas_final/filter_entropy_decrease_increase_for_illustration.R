library(tidyverse)
library(fmri.pipeline)
load("~/code/clock_analysis/coxme/fMRI_MEG_coxme_objects_no_MEDUSA_Nov23_2020")

repo_directory = "~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/f"
source(file.path(paste0(repo_directory, "get_trial_data.R")))

# get basis-wise plots

bdf <- setDT(read_csv("~/Desktop/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_trial_outputs_by_timestep.csv")) %>% 
  select(1:28) %>% pivot_longer(cols = 5:28, names_to = "time", values_to = "value") %>%
  mutate(time = parse_number(time))


trial_df <- setDT(get_trial_data(repo_directory = "~/code/clock_analysis", dataset = "mmclock_fmri"))
trial_df <- trial_df %>% select(id, run, asc_trial, run_trial, v_entropy_wi_change, v_entropy_wi_change_lag, v_entropy_wi, rt_csv, score_csv, reward) %>% arrange(id, run, run_trial) %>% group_by(id, run) %>%  
  mutate(v_entropy_wi_change_lead = lead(v_entropy_wi_change),
         v_entropy_wi_change_lead2 = lead(v_entropy_wi_change, 2),
         v_entropy_wi_change_lead3 = lead(v_entropy_wi_change, 3),
         v_entropy_wi_change_lag2 = lag(v_entropy_wi_change, 2))

df <- inner_join(bdf, trial_df, by = c("id", "asc_trial"))

#pick trials
# #   increase, then decrease                                  before entropy increase                                   after entropy increase, before decrease                   
# edf <- df %>% arrange(id, run, run_trial, time) %>% filter((v_entropy_change > 2  & v_entropy_wi_change_lead < -2) | (v_entropy_wi_change_lag > 2 & v_entropy_change < -2) |
#                                                              # after increase, before decrease
#                                                              (v_entropy_wi_change_lag2 > 2 & v_entropy_wi_change_lag < -2)                                                           )

#   decrease, then increase                                  before entropy increase                                   after entropy increase, before decrease                   
edf <- df %>% arrange(id, run, run_trial, time) %>% filter((v_entropy_wi_change <  -1.75  & v_entropy_wi_change_lead > 1.75) | (v_entropy_wi_change_lag <  -1.75 & v_entropy_wi_change > 1.75) |
                                                             # after increase, before decrease
                                                             (v_entropy_wi_change_lag2 <  -1.75 & v_entropy_wi_change_lag > 1.75)) %>%
  group_by(id, run) %>% filter(n_distinct(run_trial) >= 3)




ggplot(edf, aes(time, value, color = as.factor(run_trial), group = run_trial)) + geom_line() +
  facet_grid(id ~ run)

setwd(file.path(paste0(repo_directory, "betas_final")))

saveRDS(edf,  file = "entropy_wi_decrease_1.75z_then_increase_1.75z.rds")

ids <- unique(edf$id)
i = 2
ggplot(edf %>% filter(id == ids[i]), aes(time, value, color = as.factor(run_trial), group = run_trial)) + geom_line() +
  facet_wrap(run ~ run_trial)


ggplot(edf %>% filter(id == ids[i]), aes(run_trial, v_entropy_change)) + geom_line() +
  facet_wrap(id ~ run)

value <- c(0, 0, 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 3, 1, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0)

trial <- c(rep(1,16), rep(2, 16))
d <- as.data.frame(cbind(value, trial)) %>% group_by(trial) %>% mutate(entropy = DescTools::Entropy(value/(sum(value))))
d$time <- rep(1:16,2)       
ggplot(d, aes(time, value)) + facet_wrap(~trial) + geom_line() + geom_text(aes(x = 8, y = 2,label = entropy))
