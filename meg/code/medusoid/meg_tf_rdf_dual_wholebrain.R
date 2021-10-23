# plots *RT prediction* MLM coefficients by sensor across timepoints and frequencies
# first rsync meg_medusa_mixed_by_wholebrain.R output from the cluster

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
source("~/code/Rhelpers/theme_black.R")

repo_directory <- "~/code/clock_analysis"
data_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/rt/output"
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/rt"

clock_epoch_label = "Time relative to clock onset, seconds"
rt_epoch_label = "Time relative to outcome, seconds"
encode = T
rt_predict = F
p_adjust_method = "bonferroni"
models = c("random_slope") # random RT intercepts vs. random slopes of RT_lag and RTVmax
# models = c("random_intercept", "random_slope") # random RT intercepts vs. random slopes of RT_lag and RTVmax
print_filenames = T # print names of output files being imported
fixed_only = F # save only fixed effects statistics
reprocess = F # T = re-import output files, F = used cached .Rds objects
plots = F
by_sensor = T # make plots by sensor (otherwise, only average z stat)
noclock = T # process RT-aligned models only
freq_threshold = 40 # set to 40 for full-spectrum output
coefficients_only = F # only return df of REML coefficients, not emtrends
setwd(data_dir)
# plots ----
for (model in models) 
  message(paste0("Processing ", model))
epoch_label = "Time relative to feedback, seconds"
if (reprocess) {
  if (!noclock) {
    epoch_label = "Time relative to clock onset, seconds"
    # get clock-aligned data
    if (model=="random_intercept") {
      file_pattern <- ".*rdf.*ri.*clock"} else if (model=="random_slope") {
        file_pattern <- ".*rdf.*rs.*clock"}
    files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
    message(paste0("Found ", length(files), " files."))
    cl <- lapply(files, function(x) {
      if (print_filenames) { print(x) }
      df <- readRDS(x) 
      # Brief version, just the z-stats for plotting
      df1 <- df$emtrends_list$emt_1 %>% 
        setNames(make.names(names(.), unique = TRUE)) %>%
        select(-matches("*\\.[1-9]+$")) %>% rename(emtrend = rt_lag_sc.trend)
      df1$regressor <- "RT_t"
      df2 <- df$emtrends_list$emt_2 %>%
        setNames(make.names(names(.), unique = TRUE)) %>%
        select(-matches("*\\.[1-9]+$")) %>% rename(emtrend  = rt_vmax_lag_sc.trend)
      df2$regressor <- "RT_Vmax_t"
      df <- rbind(df1, df2) %>% dplyr::select(Freq, Time, Sensor, Pow, emtrend, std.error, reward_lag, regressor) 
      df <- df %>%
        pivot_wider(values_from=c(emtrend, std.error), names_from=c(Pow), id_cols=c(Freq, Time, Sensor, reward_lag, regressor)) %>%
        mutate(zhigh = emtrend_2/std.error_2, zlow=`emtrend_-2`/`std.error_-2`, zdiff=zhigh - zlow)
      # }
      #        df <- df %>% filter(effect=="fixed")
      return(df)
    })
    crdf <- data.table::rbindlist(cl) %>% unique() %>% rename(reward_t = reward_lag) %>%
      mutate(t = Time)#%>% distinct(Time, Freq, term, effect, group, level, .keep_all = TRUE)
    crdf$alignment <- "clock"
    message("Processed clock-aligned. \n")}
  # get RT-aligned
  if (model=="random_intercept") {
    file_pattern <- ".*rdf.*ri.*RT"} else if (model=="random_slope") {
      file_pattern <- ".*rdf.*rs.*RT"}
  # file_pattern <- "ddf_combined_entropy_rsRT|ddf_combined_entropy_change_rs_RT"
  # file_pattern <- "meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_RT|meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_finishRT"
  # file_pattern <- "entropy_rs_singleRT"
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  message(paste0("Found ", length(files), " files."))
  rl <- lapply(files, function(x) {
    if (print_filenames) { print(x) }
    df <- readRDS(x) 
    # Brief version, just the z-stats for plotting
    if (coefficients_only) {
      df <- df$coef_df_reml %>% 
        setNames(make.names(names(.), unique = TRUE)) %>% filter(effect=="fixed")
      # select(-matches("*\\.[1-9]+$")) %>%
    } else {
      df1 <- df$emtrends_list$emt_1 %>% 
        setNames(make.names(names(.), unique = TRUE)) %>%
        # select(-matches("*\\.[1-9]+$")) %>%
        rename(emtrend = rt_csv_sc.trend, reward = outcome...2)
      df1$regressor <- "RT_t"
      df2 <- df$emtrends_list$emt_2 %>%
        setNames(make.names(names(.), unique = TRUE)) %>%
        # select(-matches("*\\.[1-9]+$")) %>% 
        rename(emtrend  = rt_vmax.trend, reward = outcome...2)
      df2$regressor <- "RT_Vmax_t"
      df <- rbind(df1, df2) %>% dplyr::select(Freq, Time, Sensor, Pow, emtrend, std.error, reward, regressor)  %>%
        pivot_wider(values_from=c(emtrend, std.error), names_from=c(Pow), id_cols=c(Freq, Time, Sensor, reward, regressor)) %>%
        mutate(zhigh = emtrend_2/std.error_2, zlow=`emtrend_-2`/`std.error_-2`, zdiff=zhigh - zlow)
    }
    # }
    #        df <- df %>% filter(effect=="fixed")
    return(df)
  })
  rrdf <- data.table::rbindlist(rl) 
  # rddf$node <- sub("_group.*", "", rddf$.filename)
  rrdf$alignment <- "RT"
  if (!noclock) {offset = 4.3} else {offset = 0.3}
  rrdf <- rrdf %>% mutate(t  = Time - offset, 
                          alignment = "RT"
  ) # %>% rename(reward_t = reward)
  # saveRDS(rddf, file = "meg_ddf_wholebrain_ec_rs_rt.rds")
  message("Processed RT-aligned, merging.  \n")
  if (!noclock) {rdf <- rbind(crdf, rrdf)} else {
    rdf <- rrdf  
  }
  # deal with different frequency variable formats
  rdf$Freq <- gsub("f_", "", rdf$Freq)
  rdf$Freq <- ordered(as.numeric(substr(as.character(rdf$Freq), 1,4)))
  # deal with different sensor labels: remove leading 0s
  rdf$Sensor <- as.character(as.integer(rdf$Sensor))
  setwd(plot_dir)
  if (coefficients_only) {
    saveRDS(rdf, file = paste0("meg_rdf_wholebrain_fixef", model, ".rds"))
  } else {
    saveRDS(rdf, file = paste0("meg_rdf_wholebrain_zstats_", model, ".rds")) }
}
if (!reprocess) {
  setwd(plot_dir)
  rdf <- readRDS(paste0("meg_rdf_wholebrain_zstats_", model, ".rds"))}
# 
# ddf <- ddf  %>% ungroup() %>% 
#   filter((Time < 1.5 & Time > -1 & alignment=="rt") | (Time < 1.5 & Time > -2 & alignment=="clock")) %>% 
#   group_by(term, alignment) %>% mutate(p_fdr = p.adjust(p.value, method = p_adjust_method),
#                                        p_level_fdr = as.factor(case_when(
#                                          # p_fdr > .1 ~ '0',
#                                          # p_fdr < .1 & p_fdr > .05 ~ '1',
#                                          p_fdr > .05 ~ '1',
#                                          p_fdr < .05 & p_fdr > .01 ~ '2',
#                                          p_fdr < .01 & p_fdr > .001 ~ '3',
#                                          p_fdr <.001 ~ '4'))) %>% ungroup()
# ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
# ddf$`p, FDR-corrected` = ddf$p_level_fdr

if (plots) {
  message(paste0("Plotting ", model, " models.  \n"))
  if (!noclock) {offset = 4.3} else {offset = 0.3}
  for (outcome in c("Omission", "Reward")) {
    for (reg in c("RT_t", "RT_Vmax_t")) {
      edf <- rdf %>% filter(regressor == reg & reward_t == outcome)
      # Sensor-wise detailed plot
      if (by_sensor) {
        if (reg=="RT_t") {
          filename = (paste0("meg_tf_rt_predict_by_sensor", reg, "_", model, "_", outcome, ".pdf"))
        } else {
          filename = (paste0("meg_tf_rt_predict_by_sensor", reg, "_", model, ".pdf"))}
        pdf(file = filename, height = 20, width = 30)
        print(ggplot(edf, aes(t, Freq)) + geom_tile(aes(fill = zdiff, alpha = abs(zdiff)>2), size = .01) +
                geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
                geom_vline(xintercept = -offset + 0.3, lty = "dashed", color = "white", size = 2) +
                geom_vline(xintercept = -offset, lty = "dashed", color = "white", size = 1) +
                geom_vline(xintercept = -2, lty = "dotted", color = "grey", size = 1) +
                scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
                geom_text(data = edf, x = -offset-.1, y = 5,aes(label = "Response(t)"), size = 1, color = "white", angle = 90) +
                geom_text(data = edf, x = -offset+.4, y = 5,aes(label = "Outcome(t)"), size = 1, color = "white", angle = 90) +
                geom_text(data = edf, x = 0.5, y = 6 ,aes(label = "Clock onset (t+1)"), size = 1, color = "black", angle = 90) +
                # scale_x_continuous(limits = c(-1,1.3), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
                theme_dark() + facet_wrap(~Sensor))
        dev.off() 
      }
      # Fixed-effect (average z stat) plots for clarity
      sdf <- edf %>% select(Freq, t, zdiff) %>% group_by(Freq, t) %>% summarise(z_diff_mean = mean(zdiff)) %>% ungroup()
      if (reg=="RT_t") {
        filename = (paste0("meg_tf_rt_predict_all_", reg, "_", model, "_", outcome, ".pdf"))
      } else {
        filename = (paste0("meg_tf_rt_predict_all_", reg, "_", model, ".pdf"))
      }
      pdf(file = filename, height = 3, width = 5)
      print(ggplot(sdf, aes(t, Freq)) + geom_tile(aes(fill = z_diff_mean, alpha = abs(z_diff_mean)>2), size = .01) +
              geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
              geom_vline(xintercept = -offset + 0.3, lty = "dashed", color = "white", size = 2) +
              geom_vline(xintercept = -offset, lty = "dashed", color = "white", size = 1) +
              geom_vline(xintercept = -2, lty = "dotted", color = "grey", size = 1) +
              scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
              geom_text(data = edf, x = -offset-.1, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
              geom_text(data = edf, x = -offset+.4, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
              geom_text(data = edf, x = 0.5, y = 6 ,aes(label = "Clock onset (t+1)"), size = 2.5, color = "black", angle = 90) +
              # scale_x_continuous(limits = c(-1,1.3), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
              theme_dark() 
      )
      dev.off()
    }
  }
  
}
