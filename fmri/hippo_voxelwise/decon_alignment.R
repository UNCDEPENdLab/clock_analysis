library(tidyverse)
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/')

setwd('~/Box Sync/SCEPTIC_fMRI/long_axis_l_2.3mm/deconvolved/')
#d <- read_csv("10637_run4_long_axis_l_2.3mm.nii.gz_deconvolved.csv")
d <- read_csv ("11366_run3_long_axis_l_2.3mm_deconvolved.csv.gz")
str(d)

install.packages("oro.nifti")
mask <- oro.nifti::readNIfTI("~/code/clock_analysis/fmri/hippo_voxelwise/long_axis_l_2.3mm.nii.gz", reorient=FALSE)

mi <- which(mask > 0, arr.ind = TRUE) 
remap <- Hmisc::cut2(mask[mi], g = 20)
mask_remap <- mask
mask_remap[mi] <- as.numeric(remap)
writeNIfTI(mask_remap, "l_hippo_20bin")


trial_df <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       v_entropy_wi = scale(v_entropy),
                                       run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject

#prototype <- trial_df %>% filter(id==10637 & run==4) %>% select(run_trial, trial, clock_onset, feedback_onset, rt_csv, rt_vmax, score_csv) %>%
prototype <- trial_df %>% filter(id==11366 & run==3) %>% select(run_trial, trial, clock_onset, feedback_onset, rt_csv, rt_vmax, score_csv) %>%
  mutate(rt_csv=rt_csv/1000, rt_vmax=rt_vmax/10) %>% 
  mutate(rt_csv_cum=clock_onset + rt_csv, rt_vmax_cum=clock_onset + rt_vmax,
         rewom=factor(score_csv > 0, levels=c(FALSE, TRUE), labels=c("omission", "reward")))

prototype
str(d)

rt_vmax_cum <- prototype %>% pull(rt_vmax_cum)

dsplit <- split(d, d$vnum)

# lapply(dsplit, function(voxel) { 
#   browser()
#   res <- approx(x=voxel$time, y=voxel$decon, xout=sort(c(voxel$time, rt_vmax_cum)))
#   which_good <- res$x %in% rt_vmax_cum
#   data.frame(time=res$x[which_good], decon_evt=res$y[which_good])
# })
# 
# voxel1 <- dsplit[[1]]
# tlist <- list()
# for (t in 1:50) {
#   rt_onset <- prototype %>% filter(run_trial==t) %>% pull(rt_vmax_cum)
#   voxel1$time_map <- voxel1$time - rt_onset
#   out <- voxel1 %>% filter(time_map > -4.5 & time_map < 4.5) %>% select(time_map, time, decon) %>% mutate(run_trial=t) #%>% ggplot(aes(x=time_map, y=decon)) + geom_line()
#   tlist[[t]] <- out
# }
# 
# tdf <- bind_rows(tlist)
# ggplot(tdf, aes(x=time_map, y=decon)) + stat_smooth()
# ggplot(tdf, aes(x=time_map, y=decon)) + facet_wrap(~run_trial) + geom_line()


#scale up
time_before <- -3
time_after <- 3
results <- lapply(dsplit, function(voxel) { 
  tlist <- list()
  for (t in 1:50) {
    #rt_onset <- prototype %>% filter(run_trial==t) %>% pull(rt_vmax_cum)
    rt_onset <- prototype %>% filter(run_trial==t) %>% pull(feedback_onset)
    rewom <- prototype %>% filter(run_trial==t) %>% pull(rewom)
    magnitude <- prototype %>% filter(run_trial==t) %>% pull(score_csv)
    voxel$time_map <- voxel$time - rt_onset
    out <- voxel %>% filter(time_map > -4.5 & time_map < 4.5) %>% select(time_map, time, decon) %>% mutate(run_trial=t) #%>% ggplot(aes(x=time_map, y=decon)) + geom_line()
    to_interpolate <- voxel %>% filter(time_map > -4.5 & time_map < 4.5) %>% select(time_map, decon)
    
    #put everything onto same time grid
    interp <- approx(x=to_interpolate$time_map, y=to_interpolate$decon, xout=seq(time_before, time_after, by=1))
    interp <- interp %>% as.data.frame() %>% setNames(c("evt_time", "decon_interp")) %>% 
      mutate(run_trial=t, vnum=voxel$vnum[1], atlas_value=voxel$atlas_value[1], rewom=rewom, magnitude=magnitude)
    tlist[[t]] <- interp
  }
  
  tdf <- bind_rows(tlist)
  return(tdf)
})

allres <- bind_rows(results)
allres <- allres %>% mutate(atlas_bin=Hmisc::cut2(atlas_value, g = 4))
sumdf <- allres %>% group_by(atlas_bin, evt_time, rewom) %>% summarize(decon_interp_m=mean(decon_interp), se=plotrix::std.error(decon_interp)) %>% ungroup()

ggplot(sumdf, aes(x=evt_time, y=decon_interp_m, ymin=decon_interp_m-se, ymax=decon_interp_m+se, color=rewom)) + 
  facet_wrap(~atlas_bin) + geom_line() + geom_ribbon(alpha=0.1)

ggplot(allres %>% filter(vnum < 50), aes(x=factor(evt_time), y=decon_interp)) + 
  stat_summary(fun.data=mean_cl_boot) + 
#  geom_line(data=sumdf %>% filter(vnum < 50), aes(x=evt_time, y=decon_interp)) +
  facet_wrap(~vnum)

