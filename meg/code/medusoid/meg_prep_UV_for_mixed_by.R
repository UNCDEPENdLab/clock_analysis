library(readr)
library(ggplot2)
library(tidyr)
library(tidyverse)
basedir <- "~/code"
setwd(file.path(basedir, "clock_analysis/coxme"))
plotdir <- "~/OneDrive/collected_letters/papers/meg/plots/uv"
behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"


###############
# Prep MEG data
###############

# start with MEG stats produced by parse_sceptic_outputs.R
setwd(file.path(basedir, "clock_analysis/coxme"))
meg <- read_rds("sceptic_signals/mmclock_meg_sceptic_selective_detailedstats_Nov2020.rds")
megu <- read_rds("sceptic_signals/mmclock_meg_sceptic_fixeduv_detailedstats_Nov2020.rds")
mvdf <- meg$v_df %>% pivot_longer(-c(id, run, trial, rewFunc, y_chosen), names_to = "bin_string", values_to = "v") %>%
  mutate(bin = as.numeric(substr(bin_string, 3, 5))) %>% filter(bin<=y_chosen)
# # sanity check
# ggplot(mvdf, aes(bin, v, color = rewFunc)) + geom_smooth()
mudf <- megu$u_df %>% pivot_longer(-c(id, run, trial, rewFunc, y_chosen), names_to = "bin_string", values_to = "u") %>%
  mutate(bin = as.numeric(substr(bin_string, 3, 5))) %>% filter(bin<=y_chosen) %>% select(id, run, trial, u, bin)
# # sanity check
# ggplot(mudf, aes(bin, u, color = rewFunc)) + geom_smooth()
mvudf <- as_tibble(merge(mvdf, mudf)) %>% arrange(id, run, trial, bin) %>% mutate (
  ID = as.integer(substr(id, 1,5))
) %>% select(-id) %>% rename(id = ID)

# write time variable
mvudf$Time <- mvudf$bin/10
mvudf$response <- mvudf$y_chosen==mvudf$bin

################
# add behavioral data
trial_df <- readRDS(behavioral_data_file) %>% as.data.frame(lapply(trial_df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
})) %>% mutate(Subject = as.integer(id), Trial = trial, Run = run) %>% group_by(id, run) %>%
  mutate(abs_pe_lag = lag(abs_pe),
         v_entropy_wi_lag = lag(v_entropy_wi),
         rt_lag2_sc = lag(rt_lag_sc)
  ) %>% ungroup() 

# load("../fmri/keuka_brain_behavior_analyses/trial_df_and_vh_pe_clusters_u.Rdata")
msdf <- inner_join(trial_df, mvudf)
# ggplot(msdf, aes(bin, v, color = rewFunc)) + geom_smooth()

msdf$trial <- as.numeric(msdf$trial)
# msdf$timestep <- msdf$y_chosen
# get within-subject and within-trial value and uncertainty
msdf <- msdf %>% group_by(id, run) %>% mutate(value_wi = scale(v),
                                              uncertainty_wi = scale(u),
                                              value_b = mean(v),
                                              uncertainty_b = mean(u), 
                                              trial_neg_inv_sc = scale(-1/run_trial)) %>% ungroup() %>%
  group_by(ID, run, run_trial) %>% mutate(value_wi_t = scale(v),
                                          uncertainty_wi_t = scale(u),
                                          value_b_t = mean(v),
                                          uncertainty_b_t = mean(u)) %>% ungroup()

# filter out no-go zones at the edges
# censoring the entire first second seems too harsh, remove first 500 ms
mfdf <- msdf %>% filter(bin > 5 & bin <35)

# diagnostics on uncertainty and value distributions
setwd(plotdir)
pdf("meg_v_diagnostics.pdf", height = 6, width = 10)
ggplot(msdf %>% filter(bin >5 & bin <35), aes(timestep, value_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2)) 
dev.off()
pdf("meg_u_diagnostics.pdf", height = 6, width = 10)
ggplot(msdf %>% filter(bin >5 & bin <35), aes(timestep, uncertainty_wi_t, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2)) 
dev.off()


# ggplot(msdf %>% filter(bin >10 & bin <35), aes(time, uncertainty_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2))
# corr.test(fdf$value_wi, fdf$uncertainty_wi)
# corr.test(fdf$value, fdf$uncertainty)
# 


##Build 0/1 WV smiles (selection history with different buffer sizes)
# msdf <- msdf %>% select(ID, run, trial, bin, everything())

#use trial-wise data frame for getting lagged timesteps
trialwise <- msdf %>% filter(bin==1) %>% 
  group_by(ID, run) %>%
  mutate(timesteplag1=dplyr::lag(timestep, n=1, order_by=run_trial),
         timesteplag2=dplyr::lag(timestep, n=2, order_by=run_trial),
         timesteplag3=dplyr::lag(timestep, n=3, order_by=run_trial),
         timesteplag4=dplyr::lag(timestep, n=4, order_by=run_trial)) %>%
  ungroup() %>% select(ID, run, trial, timesteplag1, timesteplag2, timesteplag3, timesteplag4)

# msdf <- msdf %>% select(-timesteplag) %>% left_join(trialwise, by=c("ID", "run", "trial"))
msdf <- msdf %>% left_join(trialwise, by=c("ID", "run", "trial"))


get_wv_smile <- function(microdf, nlags=3, nbefore=0, nafter=1, spec=FALSE) {
  smile <- rep(0, nrow(microdf))
  lcols <- ifelse(isTRUE(spec), paste0("timesteplag", nlags), paste0("timesteplag", 1:nlags))
  wvprefix <- ifelse(isTRUE(spec), "wvs", "wv")
  vv <- microdf[1,lcols] #just first row (bin) of this trial is needed
  if (nrow(vv) > 0L) {
    allpos <- lapply(vv, function(x) {
      if (is.na(x)) {
        return(NULL)
      } else {
        return(seq.int(x-nbefore, x+nafter))
      }
    })
    
    topopulate <- unname(unlist(allpos))
    smile[topopulate[topopulate <= nrow(microdf)]] <- 1
  }
  microdf[[paste0(wvprefix, nlags, "b", nbefore, "a", nafter)]] <- smile
  return(microdf)
}

#get_wv_smile assumes bin-ordered data
msdf <- msdf %>% arrange(ID, run, trial, bin)

msdf$splitbasis <- with(msdf, paste(ID, run, trial, sep="."))

splitdf <- split(msdf, msdf$splitbasis)
splitdf <- lapply(splitdf, function(microdf) {
  microdf %>% get_wv_smile(n=3, nbefore=0, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=2) %>%
    get_wv_smile(n=4, nbefore=0, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=2) %>%
    get_wv_smile(n=1, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=2, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1, spec=TRUE)
})

mbb <- bind_rows(splitdf)

mfbb <- mbb %>% filter(bin >10 & bin <35)

#spot check
# mbb %>% filter(trial > 5) %>% select(bin, timestep, timesteplag1, timesteplag2, wv3b0a1)
#  table(mbb$wv3b0a1)
#  table(mbb$wv4b0a1)

## strategy for merging MEDUSA data: 1. write downsampled evt-time to coxme df, 2. merge in long MEDUSA df by id, run, trial, evt_time

# spot-check obsessively
# pdf('meg_check.pdf', height = 6, width = 6)
# ggplot(mbb, aes(run_trial, value_wi, color = rewFunc)) + geom_smooth()
# dev.off()
# 
# pdf('fmri_check.pdf', height = 6, width = 6)
# p1 <- ggplot(bb, aes(run_trial, value_wi, color = rewFunc)) + geom_smooth(method = 'gam')
# p2 <- ggplot(bb, aes(t1, value_wi, color = rewFunc)) + geom_smooth(method = 'gam')
# p3 <- ggplot(bb, aes(run_trial, uncertainty_wi, color = rewFunc)) + geom_smooth(method = 'gam')
# p4 <- ggplot(bb, aes(t1, uncertainty_wi, color = rewFunc)) + geom_smooth(method = 'gam')
# ggarrange(p1,p2, p3, p4)
# dev.off()

# merge in MEDUSA data -- I wonder if we have some missing censored evt_times

# use approx to interpolate

# save w/o MEDUSA
save(file = "fMRI_MEG_coxme_objects_no_MEDUSA_Nov23_2020.Rdata", bb, fbb, mbb, mfbb)

library(zoo)
if (medusa) {
  load ("~/Box/SCEPTIC_fMRI/dan_medusa/cache/clock_dan_medusa_for_coxme.Rdata")  
  clock_cox$ID <- clock_cox$id
  # medbb <- inner_join(clock_cox, bb)  %>% arrange(id, run, trial, t1)
  tmp <- left_join(clock_cox, bb)  %>% arrange(id, run, run_trial, evt_time) 
  # View(tmp[1:500,c(1:6, 31:36, 88,89, 106)])
  
  labels <- names(clock_cox[grepl("_R|_r|_L|_l", names(clock_cox))])
  medbb <-  tmp %>%  mutate_at(labels, funs(ifelse((online=="TRUE" & t1 != evt_time & !response), NA,.)))
  medbb <- medbb  %>% arrange(ID, run, run_trial,evt_time)  %>% group_by(ID, run, run_trial) %>%
    mutate_at(labels, list(interp = na.approx), na.rm = F)
  # # spot-check -- not too much variability, but OK
  # View(medbb[1:500,c(1:6, 31:35, 88,89,128, 106)])
  # # spot-check again
  # ggplot(medbb %>% filter(id==10637 & trial ==4 & run ==1), aes(t1, AIP1_R)) + 
  #   geom_point() + geom_line(aes(t1, AIP1_R_interp))
  medbb <- medbb %>% filter(!is.na(response)) %>% select(-labels)
  medfbb <- medbb %>% filter(bin > 10 & bin < 35)
}
save(file = "fMRI_MEG_coxme_objects_with_medusa_Nov24_2020.Rdata", medbb, medfbb)
