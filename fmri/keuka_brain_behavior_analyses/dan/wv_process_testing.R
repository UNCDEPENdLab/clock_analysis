library (tidyverse)
# bdf <- readRDS("/Users/hallquist/Data_Analysis/clock_analysis/fmri/mmy3_trial_df_selective_groupfixed.rds")
# load(file="/Users/hallquist/Data_Analysis/clock_analysis/coxme/clock_for_coxme_value_only_070518.RData")

source("/Users/hallquist/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/parse_sceptic_outputs.R")

setwd("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan")

#data from scratch
# sceptic_stats <- parse_sceptic_outputs("/Users/hallquist/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_sceptic_decay_fits",
#                             "/Users/hallquist/Data_Analysis/clock_analysis/fmri/behavior_files")
# 
# saveRDS(sceptic_stats, file="mmclock_fmri_sceptic_selective_detailedstats_Nov2020.rds")

# sceptic_stats <- parse_sceptic_outputs("/Users/hallquist/Data_Analysis/clock_analysis/fmri/data/mmclock_meg_sceptic_decay_fits/",
#                                        "/Users/hallquist/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/mmclock_meg")
# 
# saveRDS(sceptic_stats, file="mmclock_meg_sceptic_selective_detailedstats_Nov2020.rds")

#MEG UV
# sceptic_stats <- parse_sceptic_outputs("/Users/hallquist/Data_Analysis/clock_analysis/fmri/data/mmclock_meg_sceptic_fixeduv_fits",
#                                        "/Users/hallquist/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/mmclock_meg")
# 
# saveRDS(sceptic_stats, file="mmclock_meg_sceptic_fixeduv_detailedstats_Nov2020.rds")


#cached data for speed
sceptic_stats <- readRDS("mmclock_fmri_sceptic_selective_detailedstats_Nov2020.rds")

#build coxme-style wide data frame for plotting
fmri_wide <- sceptic_stats$v_df %>% pivot_longer(-c(id, run, trial, rewFunc, y_chosen), names_to = "bin_string", values_to = "v") %>%
  mutate(bin = as.numeric(substr(bin_string, 3, 5))) %>% select(-bin_string) %>% arrange(id, run, trial, bin)


trialwise <- fmri_wide %>% filter(bin==1) %>% 
  dplyr::rename(timestep=y_chosen) %>%
  group_by(id, run) %>%
  mutate(timesteplag1=dplyr::lag(timestep, n=1, order_by=trial),
         timesteplag2=dplyr::lag(timestep, n=2, order_by=trial),
         timesteplag3=dplyr::lag(timestep, n=3, order_by=trial),
         timesteplag4=dplyr::lag(timestep, n=4, order_by=trial)) %>%
  ungroup() %>% select(id, run, trial, timesteplag1, timesteplag2, timesteplag3, timesteplag4)

# msdf <- msdf %>% select(-timesteplag) %>% left_join(trialwise, by=c("ID", "run", "trial"))
fmri_wide <- fmri_wide %>% left_join(trialwise, by=c("id", "run", "trial"))


fmri_wide$splitbasis <- with(fmri_wide, paste(id, run, trial, sep="."))

splitdf <- split(fmri_wide, fmri_wide$splitbasis)
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

fmri_wide <- bind_rows(splitdf) %>% select(-splitbasis)


#qualitative smile plots
forplots <- fmri_wide %>% select(id, run, trial, rewFunc, bin, y_chosen, starts_with("wvs"))

pdf("tooth_plots.pdf", width=12, height=9)
for (thisid in unique(forplots$id)) {
  subdf <- forplots %>% filter(id==thisid) %>%
    dplyr::rename(tooth1=wvs1b1a1, tooth2=wvs2b1a1, tooth3=wvs3b1a1, tooth4=wvs4b1a1)
  for (thisrun in unique(subdf$run)) {
    rundf <- subdf %>% filter(run==thisrun) %>% pivot_longer(cols=starts_with("tooth"))
    gg <- ggplot(rundf, aes(x=bin, y=value, fill=name)) + geom_bar(stat="identity") + facet_wrap(~trial) +
      geom_vline(aes(xintercept=y_chosen)) + ggtitle(paste0("ID: ", rundf$id[1], ", rewFunc: ", rundf$rewFunc[1])) +
      scale_fill_brewer("History", palette="Dark2")
    plot(gg)
  }
}
dev.off()


####
#pull in betas for PEs and Entropy
load("../trial_df_and_vh_pe_clusters_u.Rdata")

betas <- df %>% filter(trial==1) %>% 
  mutate(id=as.character(id)) %>%
  select(id, pe_f1_cort_hipp, pe_f1_cort_hipp_resp, pe_f2_cerebell, pe_f3_str, pe_f3_str_resp, pe_PH_r, 
         pe_ips, pe_ips_resp, dan_h_resp, dan_l_sfg, dan_parietal, dan_r_sfg, dan_general_entropy_resp, general_entropy)

rm(df, mdf)

#what representation are we going for?
# multinomial mixed model:
#  - choose best (RTVmax, best tooth?)
#  - other tooth
#  - something else (new/explore)
#
# WV smile process unfolds by remembering recent options and their outcomes, swapping in new tooth if old
#  ones become depressed in value
#
# Need measures that capture depression of value, leading to a switch



trial_stats <- sceptic_stats$trial_stats %>% select(-rt_next, -score_next, -rt_csv, -rt_vba) %>%
  dplyr::rename(rt=y_chosen) %>%
  mutate(outcome = factor(score_csv > 0, levels=c(FALSE, TRUE), labels=c("omission", "reward"))) %>%
  group_by(id, run) %>%
  mutate(
    last_outcome = dplyr::lag(outcome, n=1, order_by=trial),
    rtlag1=dplyr::lag(rt, n=1, order_by=trial),
    rtlag2=dplyr::lag(rt, n=2, order_by=trial),
    rtlag3=dplyr::lag(rt, n=3, order_by=trial),
    rtlag4=dplyr::lag(rt, n=4, order_by=trial),
    vchosenlag1=dplyr::lag(v_chosen, n=1, order_by=trial),
    vchosenlag2=dplyr::lag(v_chosen, n=2, order_by=trial),
    vchosenlag3=dplyr::lag(v_chosen, n=3, order_by=trial),
    vchosenlag4=dplyr::lag(v_chosen, n=4, order_by=trial)) %>%
  ungroup() %>%
  mutate(
    trial=as.integer(trial),
    run_trial=case_when(
      trial >= 1 & trial <= 50 ~ trial,
      trial >= 51 & trial <= 100 ~ trial - 50L, #dplyr/rlang has gotten awfully picky about data types!!
      trial >= 101 & trial <= 150 ~ trial - 100L,
      trial >= 151 & trial <= 200 ~ trial - 150L,
      trial >= 201 & trial <= 250 ~ trial - 200L,
      trial >= 251 & trial <= 300 ~ trial - 250L,
      trial >= 301 & trial <= 350 ~ trial - 300L,
      trial >= 351 & trial <= 400 ~ trial - 350L,
      TRUE ~ NA_integer_),
    rt_swing = abs(rt - rtlag1)
  )

#use right join to force presence of betas
trial_stats <- trial_stats %>% right_join(betas)



#worker for tooth calculation
tooth_stat <- function(rdf, nbefore=1, nafter=1, lags=2, stat="old", rtlagprefix="rtlag") {
  #rt_lags <- na.omit(rdf[,paste0("rtlag", 1:lags)])
  #if (nrow(rt_lags) == 0L) { return(FALSE) } #early choice before any lags are present
  rt_lags <- na.omit(rdf[paste0(rtlagprefix, 1:lags)])
  
  if (length(rt_lags) == 0L) { return(NA) } #early choice before any lags are present
  
  rt_t <- rdf["rt"]
  rt_set <- c()
  for (rr in rt_lags) { rt_set <- c(rt_set, seq(rr-nbefore, rr+nafter)) }
  rt_set <- unique(rt_set)
  if (stat=="old") {
    rt_t %in% rt_set
  } else if (stat == "width") {
    length(rt_set)
  } else if (stat == "dist") {
    min(abs(rt_t - rt_lags))
  } else if (stat == "dist_rtvmax") {
    abs(rt_t - rdf["rt_vmax"])
  } else if (stat == "best_tooth") {
    #look at v_chosenlag vector -- does the rt fall within the tooth of the best
    vchosen_lags <- na.omit(rdf[paste0("vchosenlag", 1:lags)])
    
    distvbest <- rt_t - rt_lags[which.max(vchosen_lags)]
    best_tooth <- distvbest >= -1*nbefore & distvbest <= nafter
    if (!best_tooth) {
      if (rt_t %in% rt_set) {
        return("other tooth")
      } else {
        return("off tooth")
      }
    } else {
      return("best tooth")
    }
  }
  
}

#tooth_stat(trial_stats[3,], nbefore=2, nafter=2)

#buffer distance. Like a swing, but distance from nearest tooth
#new tooth, old tooth -- based on lag of toothness?

#variants
trial_stats$old_b2_a2_l2 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=2, nafter=2, lags=2, stat="old")
trial_stats$old_b2_a2_l3 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=2, nafter=2, lags=3, stat="old")
trial_stats$old_b2_a2_l4 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=2, nafter=2, lags=4, stat="old")

trial_stats$old_b1_a1_l2 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=2, stat="old")
trial_stats$old_b1_a1_l3 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=3, stat="old")
trial_stats$old_b1_a1_l4 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=4, stat="old")


trial_stats$width_b1_a1_l2 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=2, stat="width")
trial_stats$width_b1_a1_l3 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=3, stat="width")
trial_stats$width_b1_a1_l4 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=4, stat="width")


#distance from nearest tooth center
trial_stats$dist_b1_a1_l2 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=2, stat="dist")
trial_stats$dist_b1_a1_l3 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=3, stat="dist")
trial_stats$dist_b1_a1_l4 <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=4, stat="dist")

#hist(trial_stats$dist_b1_a1_l2)

#distance from rtvmax
trial_stats$dist_rtvmax <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, stat="dist_rtvmax")




library(tidyverse)
ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b2_a2_l2), color=rewFunc)) + geom_smooth()
ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b2_a2_l3), color=rewFunc)) + geom_smooth()
ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b2_a2_l4), color=rewFunc)) + geom_smooth()

#look at oldness, width of all teeth, and distance to nearest tooth center
ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b1_a1_l3), color=rewFunc)) + geom_smooth(method="loess")
ggplot(trial_stats, aes(x=run_trial, y=width_b1_a1_l3, color=rewFunc)) + geom_smooth(method="loess")
ggplot(trial_stats, aes(x=run_trial, y=width_b1_a1_l4, color=rewFunc)) + geom_smooth(method="loess")
ggplot(trial_stats %>% filter(dist_b1_a1_l3 > 0), aes(x=run_trial, y=dist_b1_a1_l3, color=rewFunc)) + geom_smooth(method="loess")

ggplot(trial_stats %>% filter(dist_b1_a1_l3 > 0), aes(x=rt_swing, y=dist_b1_a1_l3, color=rewFunc)) + geom_smooth(method="loess")

ggplot(trial_stats, aes(x=rt_swing, y=width_b1_a1_l4, color=rewFunc)) + geom_smooth(method="loess")

#rtvmax distance

hist(trial_stats$dist_rtvmax_b1_a1_l3)

ggplot(trial_stats %>% filter(dist_rtvmax_b1_a1_l3 > 0), aes(x=run_trial, y=dist_rtvmax_b1_a1_l3, color=rewFunc)) + geom_smooth(method="loess")


#geom_smooth(method="glm", method.args=list(family="binomial")) + geom_point()
cor.test(trial_stats$rt_swing, as.numeric(trial_stats$old_b1_a1_l3))

ggplot(trial_stats, aes(x=rt_swing, y=as.numeric(old_b1_a1_l4), color=rewFunc)) + geom_smooth(method="loess")

#add uniform random RT to data.frame to see how falloff of rt_swing scales with hitting vs. missing tooth 
trial_stats <- trial_stats %>% mutate(
  #rt_runif = round(runif(n=nrow(trial_stats), min = 0, max = 40))) %>% group_by(id, run) %>%
  rt_runif = round(grwalk(len=nrow(trial_stats), min = 0, max = 40, step_sd = 10))) %>% group_by(id, run) %>%
  mutate(
    rt_runif_lag1 = dplyr::lag(rt_runif, 1, order_by=trial),
    rt_runif_lag2 = dplyr::lag(rt_runif, 2, order_by=trial),
    rt_runif_lag3 = dplyr::lag(rt_runif, 3, order_by=trial),
    rt_runif_lag4 = dplyr::lag(rt_runif, 4, order_by=trial)) %>% ungroup() %>%
  mutate(rt_swing_runif = abs(rt_runif - rt_runif_lag1))

trial_stats$old_b1_a1_l4_runif <- apply(trial_stats %>% select(starts_with("rt")), 1, tooth_stat, nbefore=1, nafter=1, lags=4, stat="old", rtlagprefix="rt_runif_lag")
ggplot(trial_stats, aes(x=rt_swing, y=as.numeric(old_b1_a1_l4_runif), color=rewFunc)) + geom_smooth(method="loess")


#best tooth
trial_stats$besttooth_b2_a2_l3 <- factor(apply(trial_stats %>% select(starts_with(c("rt", "vchosen"))), 1, tooth_stat, nbefore=2, nafter=2, lags=3, stat="best_tooth"))

#need summary data...
sum_df <- trial_stats %>% filter(!is.na(besttooth_b2_a2_l3)) %>% 
  group_by(run_trial, rewFunc, besttooth_b2_a2_l3, pe_ips_resp) %>% 
  summarise(celln=n()) %>% group_by(run_trial, rewFunc) %>% mutate(tot=sum(celln), pct=celln/tot) %>%
  ungroup()

sum_df$besttooth_b2_a2_l3 <- ordered(sum_df$besttooth_b2_a2_l3, levels=c("off tooth", "other tooth", "best tooth"))
ggplot(sum_df, aes(x=as.factor(run_trial), y=pct, fill=besttooth_b2_a2_l3)) + 
  facet_grid(rewFunc~pe_ips_resp) + geom_bar(stat="identity")



#### DIVIDE PLOTS BY BETAS
ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b2_a2_l2), color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~pe_ips_resp)

ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b2_a2_l3), color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~pe_ips_resp)

ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b2_a2_l3), color=pe_ips_resp, lty=dan_general_entropy_resp)) + 
  geom_smooth(method="loess", se = FALSE) +
  facet_wrap(~rewFunc)

ggplot(trial_stats %>% filter(!is.na(last_outcome) & run_trial>3), aes(x=run_trial, y=as.numeric(old_b2_a2_l3), color=pe_ips_resp, lty=dan_general_entropy_resp)) + 
  geom_smooth(method="loess", se = TRUE) +
  facet_grid(last_outcome~rewFunc)

ggsave("dan_pe_oldnew.pdf", width=15, height=10)

ggplot(trial_stats %>% filter(!is.na(last_outcome)), aes(x=run_trial, y=as.numeric(old_b2_a2_l3), color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~dan_general_entropy_resp)


##

ggplot(trial_stats, aes(x=run_trial, y=as.numeric(old_b2_a2_l4), color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~pe_ips_resp)

ggplot(trial_stats, aes(x=run_trial, y=width_b1_a1_l3, color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~pe_ips_resp)

ggplot(trial_stats %>% filter(!is.na(last_outcome) & run_trial>3), aes(x=run_trial, y=width_b1_a1_l3, color=pe_ips_resp, lty=dan_general_entropy_resp)) + 
  geom_smooth(method="loess", se = FALSE) +
  facet_grid(last_outcome~rewFunc)

ggsave("dan_pe_width.pdf", width=15, height=10)



ggplot(trial_stats, aes(x=run_trial, y=width_b1_a1_l4, color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~pe_ips_resp)




ggplot(trial_stats %>% filter(dist_b1_a1_l3 > 0), aes(x=run_trial, y=dist_b1_a1_l3, color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~pe_ips_resp)


ggplot(trial_stats, aes(x=rt_swing, y=as.numeric(old_b1_a1_l4), color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~pe_ips_resp)

ggplot(trial_stats, aes(x=rt_swing, y=as.numeric(old_b1_a1_l4), color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~dan_general_entropy_resp)


library(lme4)
trial_stats$other_tooth <- as.numeric(trial_stats$besttooth_b2_a2_l3 == "other tooth")
trial_stats$best_tooth <- as.numeric(trial_stats$besttooth_b2_a2_l3 == "best tooth")
trial_stats$on_tooth <- as.numeric(trial_stats$besttooth_b2_a2_l3 != "off tooth")
trial_stats$trial_z <- as.vector(scale(trial_stats$trial))
trial_stats$inv_trial <- -1000/trial_stats$trial
trial_stats$inv_trial_z <- as.vector(scale(trial_stats$inv_trial))

#no evidence of 3-ways for best tooth or on tooth
mm1 <- glmer(best_tooth ~ (rewFunc + trial_z + pe_ips)^2 + (1|id/run), trial_stats, family=binomial,
             control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=2e5)))
summary(mm1)

mm1_e <- glmer(best_tooth ~ (rewFunc + trial_z + general_entropy)^2 + (1|id/run), trial_stats, family=binomial,
             control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=2e5)))
summary(mm1_e)


mm_other <- glmer(other_tooth ~ (rewFunc + trial_z + pe_ips)^2 + (1|id/run), trial_stats %>% filter(best_tooth==0), 
                  family=binomial, control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=2e5)))
summary(mm_other)


mm2 <- glmer(on_tooth ~ (rewFunc*trial_z*pe_ips)^2 + (1|id/run), trial_stats, family=binomial)
summary(mm2)

mm3 <- glmer(best_tooth ~ (rewFunc*inv_trial_z*pe_ips)^2 + (1|id/run), trial_stats, family=binomial)
summary(mm3)

mm4 <- glmer(on_tooth ~ (rewFunc*inv_trial_z*pe_ips)^2 + (1|id/run), trial_stats, family=binomial)
summary(mm4)


#gaussian random walk with reflecting boundaries
grwalk <- function(len, start=runif(1, min=0, max=40), step_sd=2, max=40, min=0)  {
  stopifnot(start > min && start < max) #not sure what we'd do otherwise...
  probs <- rep(NA, len)
  rvec <- rnorm(len, mean=0, sd=step_sd)
  probs[1] <- start
  for (i in 2:len) {
    test <- probs[i-1] + rvec[i]
    if (test > max || test < min) {
      incr <- -1*rvec[i]
    } else {
      incr <- rvec[i]
    }
    
    probs[i] <- probs[i-1] + incr
  }
  
  return(probs)
}

xx <- grwalk(100, step_sd = 8)
cor.test(xx, lag(xx, n=1))

plot(xx)

cor.test(trial_stats$rt, trial_stats$rtlag1)


summary(lmer(rt))

#not behaving
# ggplot(trial_stats, aes(x=as.factor(run_trial), fill=besttooth_b2_a2_l2)) + 
#   facet_wrap(~rewFunc) +  geom_bar(mapping = aes(x=as.factor(run_trial), y = ..count../sum(..count..)), stat = "count")

#SUPER SLOW!
# trial_stats22 <- trial_stats %>% 
#   pmap_dfr(function(...) {
#     rdf <- tibble(...)
#     rdf %>% mutate(tooth_stat = tooth_stat(rdf))
#   })
  


#####


# number bins
sdf <- sdf %>% group_by(ID, run, trial) %>% mutate(bin = 1:n(), time = bin/10) %>% ungroup() %>% 
  group_by(ID) %>% mutate(value_wi = scale(value),
                          uncertainty_wi = scale(uncertainty),
                          value_b = mean(value),
                          uncertainty_b = mean(uncertainty), 
                          trial_neg_inv_sc = scale(-1/trial)) %>% ungroup() %>% 
  mutate(rtlag_sc = scale(rtlag))

sdf <- sdf %>% select(ID, run, trial, bin, everything())

sdf %>% filter(ID==10637 & trial==2) %>% View()
