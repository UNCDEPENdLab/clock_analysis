# plot entropy timecourses for fixed and selective models
library(data.table)
library(tidyverse)
library(R.matlab)
library(ggbreak) 
library(patchwork)
session  = "fmri"
if (Sys.getenv("USER")=="alexdombrovski") {
  analysis_dir <- "~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan"
  setwd(analysis_dir)
  source("get_trial_data.R")
} else {
  analysis_dir <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan"
  setwd(analysis_dir)
  source("get_trial_data.R")
}

setwd(analysis_dir)

# get random uniform priors entropy
rand_df <- 
  readMat("~/OneDrive - University of Pittsburgh/Documents/skinner/projects_analyses/SCEPTIC/subject_fitting/entropy_analysis/unisession_decay_and_fixed_rand_priors_full.mat")
decay_df <- rand_df$decay.H %>% t() %>% as_tibble() %>% mutate(trial = 1:400) %>% pivot_longer(names_to = "subject", cols = (-trial)) %>%
  mutate(subject_num = parse_number(subject)) %>% rename(entropy_decay = value)
fixed_df <- rand_df$fixed.H %>% t() %>% as_tibble() %>% mutate(trial = 1:400) %>% pivot_longer(names_to = "subject", cols = (-trial)) %>%
  mutate(subject_num = parse_number(subject)) %>% rename(entropy_fixed = value)
rand_df <- merge(decay_df, fixed_df) %>% select(-subject)

rand_ldf <- rand_df %>% pivot_longer(names_to = "model", cols = c(entropy_decay, entropy_fixed))
pdf(paste0("entropy_timecourse_fmri_random_priors.pdf"), height = 1.75, width = 3)
ggplot(rand_ldf %>% filter(trial > 1), aes(trial, value, color = model)) + geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"))
dev.off()

if (session == "meg") {
  study = "mmclock_meg"
  run_length = 63} else if (session == "fmri") {study = "mmclock_fmri"
  run_length = 50}
df <- setDT(get_trial_data(repo_directory = "~/code/clock_analysis", dataset = study)) %>%
  # trial_df <- trial_df %>%
  group_by(id, run) %>%
  mutate(v_max_wi_lag = lag(v_max_wi, order_by = run_trial),
         id = as.integer(id)) %>%
  ungroup() %>%
  dplyr::rename(run_number = run) %>%
  dplyr::select(id, trial, run_number, run_trial, trial_neg_inv, rt_csv, rt_lag, v_entropy_wi, v_max_wi_lag,
                rt_vmax_lag, last_outcome, rewFunc, v_entropy, v_entropy_full)
ldf <- pivot_longer(df, cols = c(v_entropy, v_entropy_full), names_to = "model") %>%
  mutate(model = case_when(
    model=="v_entropy" ~ "information-compressing",
    model=="v_entropy_full" ~ "traditional"
  )) %>% rename(entropy = value)

# add random uniform priors
id_df <- df %>% select(id) %>% unique() %>% mutate(
  subject_num = as.numeric(1:length(unique(df$id)))
) %>% merge(rand_df, by = "subject_num")

rand_df <- id_df %>% merge(df, by = c("id", "trial"))


setwd("~/OneDrive/collected_letters/papers/meg/figures/fig_1_task_model/")
pdf(paste0("entropy_timecourse_", session, ".pdf"), height = 1.75, width = 1.75)
ggplot(ldf, aes(trial, entropy, color = model)) + 
  # geom_smooth(method = "loess", span = .1) +
  scale_color_grey() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 40, bs = "cs")) +
  scale_x_continuous(breaks = c(seq(from = 0, to = run_length*8, by = run_length))) + xlim(10, run_length*8) +
  geom_vline(xintercept = seq(from = run_length, to = run_length*8, by = run_length), lty = "dotted", size = .5) + 
  labs(x = "Trial across all runs", y = "Entropy, nats") +
  guides(color=F) +
  theme_minimal()
dev.off()
pdf(paste0("entropy_runwise_timecourse_", session, ".pdf"), height = 2, width = 2)
ggplot(ldf %>% filter(run_number>1), aes(run_trial, entropy, color = model)) + 
  # geom_smooth(method = "loess", span = .1) + 
  scale_color_grey() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 40, bs = "cs")) +
  labs(x = "Trial within run", y = "Entropy, nats") +
  # scale_x_continuous(breaks = c(seq(from = 0, to = run_length*8, by = run_length))) + xlim(10, run_length*8) +
  # geom_vline(xintercept = seq(from = run_length, to = run_length*8, by = run_length), lty = "dotted", size = .5) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme_minimal() + 
  # theme(legend.position = c(.6,.6), legend.background = element_rect(fill = "white", linetype = 0))
  theme(legend.position = "none") + 
  annotate("text", x = 30, y = 1.03, label = "Information-\ncompressing") +
  annotate("text", x = 15, y = 1.16, label = "Traditional", color = "grey40")
dev.off()

sdf <- ldf %>% group_by(trial, model) %>% summarise(median_entropy = median(entropy))
ggplot(sdf, aes(trial, median_entropy, color = model)) + geom_line() + scale_x_continuous(breaks = c(seq(from = 0, to = 504, by = 63))) +
  geom_vline(xintercept = seq(from = 63, to = 504, by = 63))

# plot with random uniform priors
rand_ldf <- pivot_longer(rand_df, cols = c(entropy_decay, entropy_fixed), names_to = "model") %>%
  mutate(model = case_when(
    model=="entropy_decay" ~ "information-compressing",
    model=="entropy_fixed" ~ "traditional"
  )) %>% rename(entropy = value)
setwd("~/OneDrive/collected_letters/papers/meg/figures/fig_1_task_model/")
pdf(paste0("entropy_timecourse_random_uniform.pdf"), height = 1.75, width = 1.75)
ggplot(rand_ldf, aes(trial, entropy, color = model)) + 
  geom_smooth(method = "loess", span = .05, size = 0.75) +
  scale_color_grey(start = 0, end = 0.6) +
  # geom_smooth(method = "gam", formula = y ~ s(x, k = 40, bs = "cs")) +
  scale_x_continuous(breaks = c(seq(from = 0, to = run_length*8, by = run_length))) + xlim(10, run_length*8) +
  geom_vline(xintercept = seq(from = run_length, to = run_length*8, by = run_length), lty = "dotted", size = .2) + 
  labs(x = "Trial across all runs", y = "Entropy, bits") +
  guides(color=F) +
  theme_minimal()
dev.off()
pdf(paste0("entropy_runwise_timecourse_random_uniform.pdf"), height = 2, width = 2)
ggplot(rand_ldf %>% filter(run_number>1), aes(run_trial, entropy, color = model)) + 
  geom_smooth(method = "loess", span = .03, size = 0.75) +
  scale_color_grey(start = 0, end = 0.6) +
  # facet_wrap(~model, scales = "free")
  # scale_y_break(c(3.07,4.39)) +
  # geom_smooth(aes(group = interaction(id, model)), se = F, size = .05, alpha = 0.03) +
  # geom_smooth(method = "gam", formula = y ~ s(x, k = 50, bs = "cs"), size = 0.7) +
  labs(x = "Trial within run", y = "Entropy, bits") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme_minimal() + 
  # theme(legend.position = c(.6,.6), legend.background = element_rect(fill = "white", linetype = 0))
  theme(legend.position = "none") + 
  annotate("text", x = 35, y = 3.3, label = "Information-\ncompressing") +
  annotate("text", x = 15, y = 4.2, label = "Traditional", color = "grey40") +
  coord_cartesian(ylim=c(2.67,4.35))
dev.off()

# pdf(paste0("entropy_runwise_timecourse_random_uniform_zoom.pdf"), height = 1, width = 1)
# ggplot(rand_ldf %>% filter(run_number>1 & model == "traditional"), aes(run_trial, entropy, color = model)) + 
#   geom_smooth(method = "loess",  span = 0.5, size = 0.75) +
#   scale_color_grey(start = 0, end = 0.6) +
#   facet_wrap(~model, scales = "free") +
#   # scale_y_break(c(3.07,4.39)) +
#   # geom_smooth(aes(group = interaction(id, model)), se = F, size = .05, alpha = 0.03) +
#   # geom_smooth(method = "gam", formula = y ~ s(x, k = 50, bs = "cs"), size = 0.7) +
#   labs(x = "Trial within run", y = "Entropy, bits") +
#   guides(color=guide_legend(override.aes=list(fill=NA))) +
#   theme_minimal() + 
#   # theme(legend.position = c(.6,.6), legend.background = element_rect(fill = "white", linetype = 0))
#   theme(legend.position = "none") + 
#   annotate("text", x = 35, y = 3.3, label = "Information-\ncompressing") +
#   annotate("text", x = 15, y = 4.2, label = "Traditional", color = "grey40") +
#   coord_cartesian(xlim = c(1, 50), ylim = c(4.39,4.405))
# dev.off()
# 

# # sanity check: make sure it holds at individual subject level
# ids <- unique(df$id)[1:10]
# ggplot(ldf %>% filter(id %in% ids), aes(run_trial, entropy, color = model)) + 
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs")) + facet_wrap(rewFunc~id)
