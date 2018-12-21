setwd("/gpfs/group/mnh5174/default/Michael/Clock_MEG")
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(foreach)
library(doParallel)
library(readr)
library(lme4)

f <- Sys.getenv('PBS_NODEFILE')
nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', 3)

cat("Node list allocated to this job\n")
print(nodelist)

cl <- makePSOCKcluster(nodelist, outfile='')
print(cl) ##; print(unclass(cl))

#clusterEvalQ(cl, { system("module use /gpfs/group/mnh5174/default/sw/modules; module load curl/7.54.0; module load r/3.4.0")
#  library(car) #the offending call
#})

registerDoParallel(cl)

#meg_files <- list.files(path="raw_data", pattern="\\.csv\\.gz", full.names=TRUE)
meg_files <- list.files(path="raw_data/clusters", pattern="\\.csv\\.gz", full.names=TRUE)


#remix <- read_csv("/gpfs/group/mnh5174/default/clock_analysis/meg/data/mmclock_meg_decay_mfx_trial_stats.csv.gz")
remix <- read_csv("/gpfs/group/mnh5174/default/clock_analysis/meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
attr(remix, "spec") <- NULL

gstats <- read_csv("/gpfs/group/mnh5174/default/clock_analysis/meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
attr(gstats, "spec") <- NULL

#ggplot(gstats, aes(x=beta_transformed, y=LL)) + geom_point() + stat_smooth() #10895 looks really bad

remix <- remix %>% separate(col=id, into=c("Subject", "VisitDate"), sep="_") %>%
  dplyr::filter(!(Subject=="11250" & VisitDate=="20140421")) %>% #duplicate data
  mutate(
    Subject=as.numeric(Subject),
    rt_csv=rt_csv/1000, #go to seconds (normalize for lmer)
    Rewarded=as.numeric(score_csv > 0)) %>% #rew/omission
  dplyr::rename(Pe=pe_max, Rt=rt_csv, Trial=trial, Run=run, Faces=emotion) %>%
  select(Subject, Trial, Run, Rt, Pe, Rewarded, Faces)
  
#meg_df <- read.csv("MEG2232_10.88398400730387.csv.gz")
#save(df_split, file="df_split.RData")
#load(file="df_split.RData")

#worker function for fitting lmer and returning essential statistics
lmer_df <- function(f, df, mname, REML=FALSE) {
  #require(DHARMa)
  #require(nortest)
  #disable REML by default to allow for AIC comparisons
  m <- tryCatch(lmer(f, df, REML=REML), error=function(e) { print(e); return(NULL) })
  if (!is.null(m)) {
    ret_df <- broom::tidy(car::Anova(m, type=3)) %>% mutate(group="anova") %>% bind_rows(broom::tidy(m))
    ret_df$Time <- df$Time[1]; ret_df$Freq <- df$Freq[1]; ret_df$model <- mname; ret_df$AIC <- AIC(m)
    #diag <- shapiro.test(resid(m)) #maxes at 5k observations
    #diag <- ad.test(resid(m))
    #ret_df$resid_W <- diag$statistic; ret_df$resid_p <- diag$p.value

    #DHARMa checks on residuals
    #dout <- simulateResiduals(fittedModel = m, refit = FALSE, n=500)
    #uniftest <- testUniformity(dout)
    #ret_df$simresid_KS_D=uniftest$statistic; ret_df$simresid_KS_p=uniftest$p.value     
    return(ret_df)
  } else {
    return(NULL)
  }
}

strict_left_join <- function(x, y, by = NULL, ...){
  by <- common_by(by, x, y)
  if(any(duplicated(y[by$y]))) {
    stop("Duplicate values in foreign key")
  } else left_join(x, y, by = by, ...)
}

allres <- list()
for (f in meg_files) {
  meg_df <- read_csv(f)
  attr(meg_df, "spec") <- NULL
  #meg_df <- read_csv("MEG2232_10.88398400730387.csv.gz")

  meg_df <- meg_df %>% mutate(Pow_dB = 10*log10(Pow + 101)) %>% arrange(Subject, Trial) %>% #-100 is the min
    filter(!Subject==11246) #this was a terrible subject in MR and the Pe distribution here is tiny (no learning?) -- need to investigate

  #initially getting more rows with left join -- turned out to be duplicate subject (11250)
  #missout <- anti_join(remix, meg_df, by=c("Subject", "Trial"))
  #missout <- anti_join(meg_df, remix, by=c("Subject", "Trial"))
  
  #merge updated PE data back into these data
  meg_df <- meg_df %>% select(-Pe, -Faces, -X1) %>% left_join(remix, by=c("Subject", "Trial"))
  
  meg_df$Pe_binary <- as.numeric(meg_df$Pe > 0) #binary representation of PPE
  meg_df$Pe_mid <- as.numeric(meg_df$Pe > quantile(meg_df$Pe, .25, na.rm=TRUE))
  #run lmer by time and frequency bin. each df contains one frequency bin
  #thus, we need to divide over time bins

  df_split <- split(meg_df, meg_df$Time)

  ## pdf("power_plots.pdf", width=10, height=10)
  ## for (tdf in df_split) {
  ##   g1 <- ggplot(tdf, aes(x=Pow)) + geom_histogram(bins=20) + theme_bw(base_size=18) + ggtitle(paste("Orig P, Time:", tdf$Time[1]))
  ##   g2 <- ggplot(tdf, aes(x=Pow_dB)) + geom_histogram(bins=20) + theme_bw(base_size=18) + ggtitle(paste("10*log10(P), Time:", tdf$Time[1]))
  ##   plot(plot_grid(g1, g2, nrow=1, align="hv"))
  ## }
  ## dev.off()
  
  #there are too many trials to wrap sanely
  # pdf("power_plots_by_trial.pdf", width=10, height=10)
  # for (tdf in df_split) {
  #   g1 <- ggplot(tdf, aes(x=Pow)) + geom_histogram(bins=20) + theme_bw(base_size=18) + ggtitle(paste("Orig P, Time:", tdf$Time[1])) + facet_wrap(~Trial)
  #   g2 <- ggplot(tdf, aes(x=Pow_dB)) + geom_histogram(bins=20) + theme_bw(base_size=18) + ggtitle(paste("10*log10(P), Time:", tdf$Time[1])) + facet_wrap(~Trial)
  #   plot(g1)
  #   plot(g2)
  # }
  # dev.off()

  #car is making workers die because of a call to curl (fixed with compilation of R 3.5.0)
  res <- foreach(tp=iter(df_split), .packages=c("lme4", "dplyr", "broom", "DHARMa", "car", "nortest"), .noexport=c("meg_df", "df_split")) %dopar% {
    tp <- tp %>% mutate(Pe_sqrt=sqrt(Pe)) %>%
      mutate_at(vars(Age, Trial, Pow_dB, Pe, Pe_sqrt), funs(z=as.vector(scale(.)))) %>% #overall z scoring (irrespective of subject)
      group_by(Subject) %>%
      mutate(Pe_wi=Pe - mean(Pe, na.rm=TRUE), Pe_pmean=mean(Pe, na.rm=TRUE), #within subject centering plus person means
        Pe_wi_z=as.vector(scale(Pe_wi))) %>% #w/i person z scoring of PEs to reduce b/w differences in magnitude
      ungroup() %>% mutate(Pe_pmean_c = Pe_pmean - mean(Pe_pmean, na.rm=TRUE),
        Pe_pmean_z = as.vector(scale(Pe_pmean)),
        Age.c = Age - mean(Age, na.rm=TRUE),
        Trial_rescale=(Trial - min(Trial))/100) #to get variance components to be on similar scales, need to lower variance of Trial (if using raw)
    
    #cov(tp %>% select(Age.c, Trial, Pow_dB, Pe), use="pairwise.complete.obs")
    #cov(tp %>% select(Age_z, Trial_z, Pow_dB), use="pairwise.complete.obs")

    #standardize predictors unless otherwise noted
    #m <- lmer(Pow_dB ~ 1 + Faces + Age + Trial + (1 + Trial | Subject), tp)    
    #m1 <- lmer(Pow_dB ~ 1 + Faces + Age_z + Trial_z + (1 + Trial_z | Subject), tp)

    #m1 <- lmer_df(Pow_dB ~ 1 + Faces + Age_z + Trial_z + (1 + Trial_z | Subject), tp, "m01")
    #m2 <- lmer_df(Pow_dB ~ 1 + Faces + Age_z * Trial_z + (1 + Trial_z | Subject), tp, "m02")
    m3 <- lmer_df(Pow_dB ~ 1 + Faces * Age_z + Trial_z + (1 + Trial_z | Subject), tp, "baseline")
    m3_run <- lmer_df(Pow_dB ~ 1 + Faces * Age_z + Trial_z + (1 + Trial_z | Subject/Run), tp, "baseline_run")
    #m4 <- lmer_df(Pow_dB ~ 1 + Faces * Age_z * Trial_z + (1 + Trial_z | Subject), tp, "m04")

    #str(anova(m1, m2, m3, m4))

    robust_baseline <- lmer_df(Pow_dB ~ 1 + Faces * Age_z + Trial_z + Rewarded + (1 + Trial_z | Subject/Run), tp, "robust_baseline")
    
    #add PE in (models 1-4 do not have it)
    pe_baseline <- lmer_df(Pow_dB ~ 1 + Faces + Pe_z + Trial_z + (1 | Subject), tp, "pe_baseline")
    pe_baseline_rewarded <- lmer_df(Pow_dB ~ 1 + Faces + Pe_z + Trial_z + Rewarded + (1 | Subject), tp, "pe_baseline_rewarded")
    pe_binary <- lmer_df(Pow_dB ~ 1 + Faces + Pe_binary + Trial_z + (1 | Subject), tp, "pe_binary")
    #pe_q25 <- lmer_df(Pow_dB ~ 1 + Faces + Pe_mid + Trial_z + (1 | Subject), tp, "pe_q25")

    #per Alex: would be nice to look at (1|Subject/Run) [nesting] Don't have run at the moment
    
    m5 <- lmer_df(Pow_dB ~ 1 + Faces * Pe_z + Faces * Age_z + Trial_z + (1 + Trial_z | Subject), tp, "m05")
    m6 <- lmer_df(Pow_dB ~ 1 + Pe_z * Age_z * Faces + Trial_z + (1 + Trial_z | Subject), tp, "m06")
    m7 <- lmer_df(Pow_dB ~ 1 + Pe_z * Age_z * Faces + Trial_z + (1 + Trial_z + Pe_z | Subject), tp, "m07")
    m8 <- lmer_df(Pow_dB ~ 1 + Pe_z * Age_z * Faces * Trial_z + (1 + Trial_z + Pe_z | Subject), tp, "m08")
    #m9 <- lmer_df(Pow_dB ~ 1 + Pe_sqrt_z * Age_z * Faces * Trial_z + (1 + Trial_z + Pe_sqrt_z | Subject), tp, "m09") #same as m8, but with sqrt Pe
    
    m10 <- lmer_df(Pow_dB ~ 1 + Pe_wi * Age_z * Faces + Pe_pmean_c * Age_z * Faces + Trial_z + (1 + Trial_z | Subject), tp, "m10") #add wi versus between as fixed
    m11 <- lmer_df(Pow_dB ~ 1 + Pe_wi_z * Age_z * Faces + Pe_pmean_z * Age_z * Faces + Trial_z + (1 + Trial_z | Subject), tp, "m11") #add z-scored wi versus between as fixed
    #m12 <- lmer_df(Pow_dB ~ 1 + Trial_z + Pe_z + Age_z + Age_z*Pe_z + Faces + Faces*Pe_z + Faces*Age_z*Pe_z + (1 + Pe_z | Subject), tp, "m12") #based on Kai's model (just z scored and using dB)

    #currently problematic due to tiny variance component for pe_wi
    #m10 <- lmer_df(Pow_dB_z ~ 1 + Pe_wi * Age_z * Faces + Pe_pmean_c * Age_z * Faces + Trial_z + (1 + Trial_z + Pe_wi | Subject), tp, "m10")

    #use z-scored variables for completely standardized solution (gets rid of convergence complaints)
    #m <- lmer(Pow_dB_z ~ 1 + Faces + Age_z + Trial_z + (1 + Trial_z | Subject), tp)
    #m2 <- lmer(Pow_dB_z ~ 1 + Faces * Age_z + Trial_z + (1 + Trial_z | Subject), tp)
    #m2 <- lmer(Pow_dB_z ~ 1 + Faces * Age_z + Trial_z + Pe + (1 + Trial_z | Subject), tp)
    
    #rbind(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12)
    #rbind(m3, pe_baseline, pe_binary, pe_q25, m5, m6, m7, m8, m9, m10, m11)
    rbind(m3, m3_run, robust_baseline, pe_baseline, pe_baseline_rewarded, pe_binary, m5, m6, m7, m8, m10, m11)
  }

  #allres[[ as.character(unique(meg_df$Freq)[1]) ]] <- res
  #allres[[ make.names(f) ]] <- res
  save(file=file.path("output", gsub(".csv.gz", "_lmerfits.RData", basename(f), fixed=TRUE)), res)
  #sink(paste0("warnings_", meg_df$Freq[1], ".txt"))
  #print(warnings())
  #sink()
}

#save(allres, file="lmer_results_m1-m9_clusters.RData")
stopCluster(cl)

