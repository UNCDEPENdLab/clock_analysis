##visualize hippo decon time series
library(tidyverse)
library(forecast)
setwd(file.path(getMainDir(), "clock_analysis/fmri/hippo_voxelwise/deconv_visualize"))
trial_df <- read_csv(file.path(getMainDir(), "clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz"))
trial_df <- trial_df %>%
  group_by(id, run) %>%  
  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
                rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                rt_lag = lag(rt_csv) ,
                omission_lag = lag(score_csv==0),
                rt_vmax_lag = lag(rt_vmax),
                v_entropy_wi = scale(v_entropy),
                run_trial=case_when(
                  trial >= 1 & trial <= 50 ~ trial,
                  trial >= 51 & trial <= 100 ~ trial - 50, #dplyr/rlang has gotten awfully picky about data types!!
                  trial >= 101 & trial <= 150 ~ trial - 100,
                  trial >= 151 & trial <= 200 ~ trial - 150,
                  trial >= 201 & trial <= 250 ~ trial - 200,
                  trial >= 251 & trial <= 300 ~ trial - 250,
                  trial >= 301 & trial <= 350 ~ trial - 300,
                  trial >= 351 & trial <= 400 ~ trial - 350,
                  TRUE ~ NA_real_)) %>% ungroup() %>%
  dplyr::mutate(rt_csv=rt_csv/1000, rt_vmax=rt_vmax/10) %>% 
  mutate(rt_vmax_cum=clock_onset + rt_vmax)

decon_files  <- list.files(path="data", pattern="*.csv.gz", full.names=TRUE)

dlist <- list()
for (d in 1:length(decon_files)) {
  
  df <- read.csv(decon_files[d], stringsAsFactors = FALSE)
  df$axis_bin <- cut(df$atlas_value, c(0, .2, .5, .75, 1))
  df <- filter(df, axis_bin %in% c("(0,0.2]", "(0.75,1]")) %>% mutate(time=time-1) %>%
    select(-x, -y, -z, -vnum) %>% group_by(atlas_name, subid, run_num, contingency, emotion, time, axis_bin) %>%
    summarize(decon=mean(decon)) %>% ungroup() %>% rename(id=subid, run=run_num, rewFunc=contingency) %>% droplevels()
  
  dlist[[d]] <- df
  
}

alldf <- dplyr::bind_rows(dlist)

combined_df <- trial_df %>% select(id, run, trial, rewFunc, emotion, clock_onset, feedback_onset, iti_onset, v_entropy, pe_max, rt_csv) %>%
  inner_join(alldf)

#ggplot(alldf, aes(x=time, y=decon, color=atlas_name, lty=axis_bin)) + geom_line() + facet_grid(id~run)
#ggplot(alldf %>% filter(atlas_name=="long_axis_l_2.3mm.nii.gz"), aes(x=time, y=decon, color=axis_bin)) + geom_line() + facet_grid(id~run)

alldf <- alldf %>% mutate(side=if_else(atlas_name=="long_axis_l_2.3mm.nii.gz", "l", "r"), pa=if_else(axis_bin=="(0,0.2]", "post", "ant"), side_pa=paste(side, pa, sep="_"))

trialdf_subset <- trial_df %>% filter(id %in% unique(alldf$id))
dfsplit <- split(alldf, alldf$id)
tdfsplit <- split(trialdf_subset, trialdf_subset$id)
names(dfsplit)==names(tdfsplit)

pdf("vis_decon.pdf", width=24, height=15)
for (this_subj in 1:length(dfsplit)) {
  events_to_plot <- tdfsplit[[this_subj]] %>% select(id, run, clock_onset, feedback_onset) %>% gather(key="event", value="evt_time", clock_onset, feedback_onset)
  
  g <- ggplot(dfsplit[[this_subj]], aes(x=time, y=decon, color=side_pa)) + geom_line() + facet_wrap(~run, ncol=1, scales="free_x") +
    geom_vline(data = events_to_plot, aes(lty=event, xintercept=evt_time)) + ggtitle(paste0("Subject: ", names(dfsplit)[this_subj])) +
    theme_bw(base_size=14)
  
  plot(g)
}
dev.off()


pdf("ccf_decon.pdf", width=20, height=12)
dfsplit <- split(alldf, alldf$id)
for (this_subj in 1:length(dfsplit)) {
  par(mfrow=c(4,2))
  df_wide <- dfsplit[[this_subj]] %>% droplevels() %>% na.omit() %>% select(-axis_bin) %>% spread(key=pa, value=decon) %>% ungroup()
  df_wide_split <- split(df_wide, df_wide$run)
  for (toplot in df_wide_split) {
    #lplot <- toplot %>% filter(side==)
    #g1 <- ggCcf(, y, lag.max = NULL, type = c("correlation"), plot = TRUE)
    browser()
  }
  ccf(df_wide_split[[1]]$`(0,0.2]`, df_wide_split[[1]]$`(0.75,1]`)
  p <- 
  g <- ggplot(dfsplit[[this_subj]] %>% filter(atlas_name=="long_axis_l_2.3mm.nii.gz"), aes(x=time, y=decon, color=axis_bin)) + geom_line() + facet_wrap(~run, ncol=1, scales="free_x") +
    geom_vline(data = events_to_plot, aes(lty=event, xintercept=evt_time))
  
  plot(g)
}
dev.off()