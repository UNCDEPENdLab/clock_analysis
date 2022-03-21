# inspect combined effect for reward theta desync and echange late beta desync

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
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/"

reward <- readRDS("/Users/alexdombrovski/Downloads/meg_reward_sensor_subject_ranefs.Rds")
# lb = late beta
lb <- readRDS("/Users/alexdombrovski/Downloads/meg_echange_sensor_subject_ranefs (1).Rds")
lb <- lb %>% filter(Time > .5 & Time < .8 & (Freq == 11.8 | Freq == 14.1 | Freq == 16.8))
# inspect entropy change
setwd(paste0(plot_dir, "entropy_change"))


pdf("echange_sensor_sub_ran_coefs.pdf", width = 8, height = 6)
ggplot(lb, aes(subject_ran_coefs, sensor_ran_coefs)) + geom_point()
dev.off()

pdf("echange_sensor_sub_ran_vals.pdf", width = 8, height = 6)
ggplot(lb, aes(subject_ran_vals, sensor_ran_vals)) + geom_point()
dev.off()

h1 <- ggplot(lb, aes(subject_ran_vals)) + geom_histogram()
h2 <- ggplot(lb, aes(sensor_ran_vals)) + geom_histogram()
h3 <- ggplot(lb, aes(fixed_effect)) + geom_histogram()
h4 <- ggplot(lb, aes(combined_effect)) + geom_histogram()

pdf("echange_sensor_sub_ran_coef_histograms.pdf", width = 6, height = 10)
print(ggarrange(h1, h2, h3, h4, ncol = 1))
dev.off()


# rth = reward theta
setwd(paste0(plot_dir, "reward"))
rth <- reward %>% filter(Time > .2 & Time < .4 & (Freq == 5 | Freq == 5.9 | Freq == 7))

pdf("reward_sensor_sub_ran_vals.pdf", width = 8, height = 6)
ggplot(lb, aes(subject_ran_vals, sensor_ran_vals)) + geom_point()
dev.off()

r1 <- ggplot(rth, aes(subject_ran_vals)) + geom_histogram()
r2 <- ggplot(rth, aes(sensor_ran_vals)) + geom_histogram()
r3 <- ggplot(rth, aes(fixed_effect)) + geom_histogram()
r4 <- ggplot(rth, aes(combined_effect)) + geom_histogram()
pdf("reward_sensor_sub_ran_coef_histograms.pdf", width = 6, height = 10)
print(ggarrange(r1, r2, r3, r4, ncol = 1))
dev.off()


