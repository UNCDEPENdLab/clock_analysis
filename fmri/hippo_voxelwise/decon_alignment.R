library(tidyverse)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise"))

#function to get interpolated event locked data
event_lock_decon <- function(d_by_bin, trial_df, event="feedback_onset", time_before=-3, time_after=3) {
  results <- lapply(d_by_bin, function(bin_data) {
    subj_df <- trial_df %>% filter(id==bin_data$subid[1] & run==bin_data$run_num[1])
    trials <- subj_df %>% pull(run_trial) %>% sort()

    tlist <- list()
    for (t in trials) {
      evt_onset <- subj_df %>% filter(run_trial==t) %>% pull(!!event)
      
      if (is.na(evt_onset)) { next } #times are missing for some events on early trials (like rt_vmax_cum)
      bin_data$time_map <- bin_data$time - evt_onset
      
      #add time on either side of interpolation grid to have preceding time points that inform linear interp
      to_interpolate <- bin_data %>% filter(time_map > time_before - 1.5 & time_map < time_after + 1.5) %>% 
        select(time_map, decon, vnum) %>% group_by(time_map) %>%
        summarize(vox_sd=sd(decon, na.rm=TRUE), decon=mean(decon, na.rm=TRUE))
      
      #for checking heterogeneity
      #ggplot(to_interpolate, aes(x=time_map, y=decon, color=factor(vnum))) + geom_line()

      #rare, but if we have no data at tail end of run, we may not be able to interpolate
      if (nrow(to_interpolate) < 2) {
        cat("For subject: ", subj_df$id[1], ", insufficient interpolation data for run: ", subj_df$run[1], ", trial: ", t, "\n", file="evtlockerrors.txt", append=TRUE, sep="")
        next
      }
      
      #put everything onto same time grid
      interp <- approx(x=to_interpolate$time_map, y=to_interpolate$decon, xout=seq(time_before, time_after, by=1)) %>%
        as.data.frame() %>% setNames(c("evt_time", "decon_interp"))
      interp_sd <- approx(x=to_interpolate$time_map, y=to_interpolate$vox_sd, xout=seq(time_before, time_after, by=1)) %>%
        as.data.frame() %>% setNames(c("evt_time", "sd_interp"))
      
      interp_all <- interp %>% left_join(interp_sd, by="evt_time") %>% 
        mutate(run_trial=t, axis_bin=bin_data$axis_bin[1])
      
      tlist[[t]] <- interp_all
    }
    
    tdf <- bind_rows(tlist)
    tdf$id <- subj_df$id[1] #tack on identifying columns
    tdf$run <- subj_df$run[1]
    return(tdf)
  })
  
  results <- bind_rows(results) %>% arrange(id, run, run_trial, evt_time) %>%
    mutate(event=!!event) %>% select(id, run, run_trial, event, evt_time, axis_bin, everything())
  
  return(results)
}


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
    

base_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise"

atlas_dirs <- list.dirs(file.path(base_dir, "deconvolved_timeseries_unsmoothed"), recursive=FALSE)

cl <- makeCluster(40)
registerDoParallel(cl)
on.exit(stopCluster(cl))

events <- c("clock_onset", "feedback_onset", "rt_vmax_cum", "clock_long", "feedback_long")
nbins <- 12 #splits along axis

for (a in atlas_dirs) {
  aname <- basename(a)

  mask <- oro.nifti::readNIfTI(file.path(getMainDir(), "clock_analysis/fmri/hippo_voxelwise", aname), reorient=FALSE)
  mi <- which(mask > 0, arr.ind = TRUE)
  
  if (length(unique(mask[mi])) == 1L) {
    bin_cuts <- seq(min(mask[mi]) - 1e-5, max(mask[mi]) + 1e-5, length.out=2) #for masks without continuous gradient
  } else {
    bin_cuts <- seq(min(mask[mi]) - 1e-5, max(mask[mi]) + 1e-5, length.out=nbins+1)
  }
  
  afiles <- list.files(file.path(a, "deconvolved"), full.names = TRUE)

  for (e in events) {
    if (e == "clock_long") {
      evt_col <- "clock_onset"
      time_before=-1
      time_after=10
    } else if (e == "feedback_long") {
      evt_col <- "feedback_onset"
      time_before=-1
      time_after=10
    } else {
      evt_col <- e
      time_before=-3
      time_after=3
    }

    out_name <- file.path("deconvolved_evt_locked_unsmoothed", paste0(aname, "_", e, "_decon_locked.csv.gz"))
    if (file.exists(out_name)) {
      message("Output file already exists: ", out_name)
      next
    }

    elist <- foreach(fname=iter(afiles), .packages = c("dplyr", "readr")) %dopar% {
      d <- read_csv(fname) %>% mutate(axis_bin=cut(atlas_value, bin_cuts)) %>%
        select(-atlas_name, -x, -y, -z)
      dsplit <- split(d, d$axis_bin)

      if (all(is.na(d$decon))) { browser() }

      subj_lock <- tryCatch(event_lock_decon(dsplit, trial_df, event = evt_col, time_before=time_before, time_after=time_after),
        error=function(err) { cat("Problems with event locking ", fname, " for event: ", e, "\n  ", as.character(err), "\n\n", file="evtlockerrors.txt", append=TRUE); return(NULL) })
      return(subj_lock)
    }
    
    all_e <- bind_rows(elist)
    all_e$atlas <- aname
    write_csv(all_e, path=out_name)
  }
}

#Times in the deconvolved files should reflect the +2 seconds for the dropped volumes
#So, this should align with the rt_csv columns appropriately.


# 
# setwd('~/Box/SCEPTIC_fMRI/long_axis_l_2.3mm/deconvolved/')
# #d <- read_csv("10637_run4_long_axis_l_2.3mm.nii.gz_deconvolved.csv")
# d <- read_csv ("11366_run3_long_axis_l_2.3mm_deconvolved.csv.gz")
# str(d)
# 
# dirfiles <- list.dirs()
# 
# mask <- oro.nifti::readNIfTI("~/code/clock_analysis/fmri/hippo_voxelwise/long_axis_l_2.3mm.nii.gz", reorient=FALSE)
# 
# # mi <- which(mask > 0, arr.ind = TRUE) 
# # remap <- Hmisc::cut2(mask[mi], g = 20)
# # mask_remap <- mask
# # mask_remap[mi] <- as.numeric(remap)
# # writeNIfTI(mask_remap, "l_hippo_20bin")
# 
# 
# 
# #prototype <- trial_df %>% filter(id==10637 & run==4) %>% select(run_trial, trial, clock_onset, feedback_onset, rt_csv, rt_vmax, score_csv) %>%
# prototype <- trial_df %>% filter(id==11366 & run==3) %>% select(run_trial, trial, clock_onset, feedback_onset, rt_csv, rt_vmax, score_csv) %>%
#   mutate(rt_csv=rt_csv/1000, rt_vmax=rt_vmax/10) %>% 
#   mutate(rt_csv_cum=clock_onset + rt_csv, rt_vmax_cum=clock_onset + rt_vmax,
#          rewom=factor(score_csv > 0, levels=c(FALSE, TRUE), labels=c("omission", "reward")))
# 
# # prototype
# # str(d)
# 
# rt_vmax_cum <- prototype %>% pull(rt_vmax_cum)
# 
# dsplit <- split(d, d$vnum)
# 
# # lapply(dsplit, function(voxel) { 
# #   browser()
# #   res <- approx(x=voxel$time, y=voxel$decon, xout=sort(c(voxel$time, rt_vmax_cum)))
# #   which_good <- res$x %in% rt_vmax_cum
# #   data.frame(time=res$x[which_good], decon_evt=res$y[which_good])
# # })
# # 
# # voxel1 <- dsplit[[1]]
# # tlist <- list()
# # for (t in 1:50) {
# #   rt_onset <- prototype %>% filter(run_trial==t) %>% pull(rt_vmax_cum)
# #   voxel1$time_map <- voxel1$time - rt_onset
# #   out <- voxel1 %>% filter(time_map > -4.5 & time_map < 4.5) %>% select(time_map, time, decon) %>% mutate(run_trial=t) #%>% ggplot(aes(x=time_map, y=decon)) + geom_line()
# #   tlist[[t]] <- out
# # }
# # 
# # tdf <- bind_rows(tlist)
# # ggplot(tdf, aes(x=time_map, y=decon)) + stat_smooth()
# # ggplot(tdf, aes(x=time_map, y=decon)) + facet_wrap(~run_trial) + geom_line()
# 
# 
# #scale up
# event_lock_decon <- function(d_by_split, trial_df, event="feedback_onset", time_before=-3, time_after=3) {
#   
# }
# 
# results <- lapply(dsplit, function(voxel) { 
#   tlist <- list()
#   for (t in 1:50) {
#     #rt_onset <- prototype %>% filter(run_trial==t) %>% pull(rt_vmax_cum)
#     rt_onset <- prototype %>% filter(run_trial==t) %>% pull(feedback_onset)
#     rewom <- prototype %>% filter(run_trial==t) %>% pull(rewom)
#     magnitude <- prototype %>% filter(run_trial==t) %>% pull(score_csv)
#     voxel$time_map <- voxel$time - rt_onset
#     out <- voxel %>% filter(time_map > -4.5 & time_map < 4.5) %>% select(time_map, time, decon) %>% mutate(run_trial=t) #%>% ggplot(aes(x=time_map, y=decon)) + geom_line()
#     to_interpolate <- voxel %>% filter(time_map > -4.5 & time_map < 4.5) %>% select(time_map, decon)
#     
#     #put everything onto same time grid
#     interp <- approx(x=to_interpolate$time_map, y=to_interpolate$decon, xout=seq(time_before, time_after, by=1))
#     interp <- interp %>% as.data.frame() %>% setNames(c("evt_time", "decon_interp")) %>% 
#       mutate(run_trial=t, vnum=voxel$vnum[1], atlas_value=voxel$atlas_value[1], rewom=rewom, magnitude=magnitude)
#     tlist[[t]] <- interp
#   }
#   
#   tdf <- bind_rows(tlist)
#   return(tdf)
# })
# 
# allres <- bind_rows(results)
# allres <- allres %>% mutate(atlas_bin=Hmisc::cut2(atlas_value, g = 4))
# sumdf <- allres %>% group_by(atlas_bin, evt_time, rewom) %>% summarize(decon_interp_m=mean(decon_interp), se=plotrix::std.error(decon_interp)) %>% ungroup()
# 
# ggplot(sumdf, aes(x=evt_time, y=decon_interp_m, ymin=decon_interp_m-se, ymax=decon_interp_m+se, color=rewom)) + 
#   facet_wrap(~atlas_bin) + geom_line() + geom_ribbon(alpha=0.1)
# 
# ggplot(allres %>% filter(vnum < 50), aes(x=factor(evt_time), y=decon_interp)) + 
#   stat_summary(fun.data=mean_cl_boot) + 
# #  geom_line(data=sumdf %>% filter(vnum < 50), aes(x=evt_time, y=decon_interp)) +
#   facet_wrap(~vnum)
# 
# 
# #LEFTOVERS
# 
# #xx <- cut(mask[mi], breaks=bin_cuts)
# 
# #just get cutpoints to apply elsewhere (this gets flaky)
# #remap <- Hmisc::cut2(mask[mi], g = nbins) #check
# #cuts <- Hmisc::cut2(mask[mi], g = nbins, onlycuts = TRUE)
# #xx <- cut(mask[mi], breaks=cuts)
# #sum(mask[mi] >= cuts[4] & mask[mi] < cuts[5])
# #table(remap)[4]
# 
# 
