# review and replot simulations from Cognition paper
# contrast exploitation between full and selective


library(tidyverse)
library(lme4)
# library(ggpubr)
library(car)
library(viridis)
# library(ggnewscale)
# library(RColorBrewer)
library(psych)
# library(corrplot)
source("~/code/Rhelpers/theme_black.R")

repo_directory <- "~/code/clock_analysis"
data_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/simulations/simdf"
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots"

print_filenames = T

# read in the two models, full and selective
models <- c("fixedLR_softmax", "fixedLR_decay", "null", "fixed_uv", "kalman_softmax", "kalman_uv_sum")
model_results <- lapply(models, function(mm) {
  setwd(file.path(data_dir))
  file_pattern <- file.path(paste0("simdf_.*", mm, "_.*.csv"))
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  message(paste0("Found ", length(files), " files."))
  cl <- lapply(files, function(x) {
    if (print_filenames) { print(x) }
    df <- readr::read_csv(x) %>% mutate(model = mm)
    return(df)
  })
  mdf <- data.table::rbindlist(cl)
  return(mdf)
})
mdf <- as_tibble(data.table::rbindlist(model_results)) %>% filter(trial<51 & !str_detect(model, "kalman")) %>% group_by(model, rewFunc, run) %>% arrange(trial) %>%
  mutate(rt_lag = lag(rt),
         rt_swing = abs(rt - rt_lag))
ggplot(mdf, aes(trial, rt_swing, color = model)) + geom_smooth(method = "loess")
summary_df <- mdf %>% group_by(model, trial) %>% summarise(median_rt_swing = median(rt_swing),
                                                           mean_score = mean(score))
ggplot(summary_df, aes(trial, median_rt_swing, color = model)) + geom_smooth(se = F) + theme_minimal()
ggplot(summary_df, aes(trial, mean_score, color = model)) + geom_smooth(se = F)

ggplot(mdf, aes(trial, score, color = model)) + geom_smooth(method = "loess", se = F)
m1 <- lmer(scale(score) ~ model + rewFunc + scale(trial) + (1|rewFunc), mdf)
summary(m1)
m2 <- lmer(scale(rt_swing) ~ model*scale(trial) + rewFunc + (1|rewFunc), mdf)
summary(m2)

sdf <- mdf %>% group_by(model) %>% summarise(mean_return = mean(score))
ggplot(sdf, aes(model, mean_return)) + geom_point()
