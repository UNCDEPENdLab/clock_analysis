# 2021-03-02 AndyP
# This function creates beta weight by significance plots
# get <x,y,z> coordinates of vmPFC mask regions
# 


library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(stargazer)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(car)

clock_folder <- "/Users/andypapale/clock_analysis" # Andrew


setwd('~/vmPFC/')

L <- list.files(pattern = '.Brain_to_Behavior.')
b2b_df <- data.frame(effect=character(),group=character(),term=character(),
                     estimate=double(),std.error=double(),statistic=double(),
                     namevar=character(),betavar=character(),betastat=double())
b2b_mdf <- data.frame(effect=character(),group=character(),term=character(),
                     estimate=double(),std.error=double(),statistic=double(),
                     namevar=character(),betavar=character(),betastat=double())
betatstat <- NULL
for (i in 1:length(L)){
  load(L[i])
  
  # get betavar
  # if (i==1){
  #   betavar = "pe-clock"
  # } else if (i==2){
  #   betavar = "pe-feedback"
  # } else if (i==3){
  #   betavar = "rtvmax-clock"
  # } else if (i==4){
  #   betavar = "rtvmax-feedback"
  # } else if (i==5){
  #   betavar = "vchosen-clock"
  # } else if (i==6){
  #   betavar = "vchosen-feedback"
  # } else if (i==7){
  #   betavar = "ventropy-clock"
  # } else if (i==8){
  #   betavar = "ventropy-feedback"
  # }
  
  # dfid <- unique(df$id)
  # # reconstruct betas for first model
  # id <- NULL
  # pfc11aR <- NULL
  # pfc11b47R <- NULL
  # pfc11L <- NULL
  # pfc11m13L <- NULL
  # pfc14m3225L <- NULL
  # pfc14m3225R <- NULL
  # pfc14r14c11mL <- NULL
  # pfc14r14c11mR <- NULL
  # pfc2432R <- NULL
  # pfc47L <- NULL
  # pfcd10R <- NULL
  # pfcfp10L <- NULL
  # pfcfp10R <- NULL
  # pfcp3224L <- NULL
  # pfcrdb10L <- NULL
  # pfcrl10R <- NULL
  # for (j in 1:length(dfid)){
  #   idix <- which(df$id==dfid[j])
  #   ix <- idix[1]
  #   id[j] <- dfid[j]
  #   pfc11aR[j] <- df$pfc11aR[ix]
  #   pfc11b47R[j] <- df$pfc11b47R[ix]
  #   pfc11L[j] <- df$pfc11L[ix]
  #   pfc11m13L[j] <- df$pfc11m13L[ix]
  #   pfc14m3225L[j] <- df$pfc14m3225L[ix]
  #   pfc14m3225R[j] <- df$pfc14m3225R[ix]
  #   pfc14r14c11mL[j] <- df$pfc14r14c11mL[ix]
  #   pfc14r14c11mR[j] <- df$pfc14r14c11mR[ix]
  #   pfc2432R[j] <- df$pfc2432R[ix]
  #   pfc47L[j] <- df$pfc47L[ix]
  #   pfcd10R[j] <- df$pfcd10R[ix]
  #   pfcfp10L[j] <- df$pfcfp10L[ix]
  #   pfcfp10R[j] <- df$pfcfp10R[ix]
  #   pfcp3224L[j] <- df$pfcp3224L[ix]
  #   pfcrdb10L[j] <- df$pfcrdb10L[ix]
  #   pfcrl10R[j] <- df$pfcrl10R[ix]
  # }
  # 
  # h <- data.frame(id,pfc11aR,pfc11b47R,pfc11L,pfc11m13L,pfc14m3225L,pfc14m3225R,
  #                 pfc14r14c11mL,pfc14r14c11mR,pfc2432R,pfc47L,pfcd10R,pfcfp10L,
  #                 pfcfp10R,pfcp3224L,pfcrdb10L,pfcrl10R)
  # 
  # h<-reshape2::melt(h,id.vars=c("id")) #tiyverse pivot_wider
  # h<-h %>% dplyr::rename("region"=variable,"beta"=value)
  
  if (i==1){
    load('2021-03-16-Schaeffer-vmPFC-pe_max.Rdata')
    betavar = 'pe_max'
  } else if (i==2) {
    load('2021-03-16-Schaeffer-vmPFC-rt_vmax_change.Rdata')
    betavar = 'rt_vmax_change'
  } else if (i==3) {
    load('2021-03-16-Schaeffer-vmPFC-v_chosen.Rdata')
    betavar = 'v_chosen'
  } else if (i==4) {
    load('2021-03-16-Schaeffer-vmPFC-v_entropy.Rdata')
    betavar = 'v_entropy'
  }
  
  h <- h %>% mutate(side = substr(Schaefer_label, nchar(Schaefer_label), nchar(Schaefer_label)),
                    region = substr(Schaefer_label, 1, nchar(Schaefer_label)))
  
  
  m1 <- lmer(beta ~ region + (1|id), h)
  m1 <- tidy(m1)
  
  pfc_df <- grep(pattern='pfc',colnames(df))
  for (k in 1:length(pfc_df)){
    tempdf <- df %>% select(id,run,rt_csv_sc,trial_neg_inv_sc,rt_lag_sc,rt_vmax_lag_sc,run_trial,rt_vmax,rewFunc,
                            last_outcome,v_max_wi_lag,v_entropy_wi,pfc_df[k])
    
  if (k==1){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11m13L)
    namevar = "pfc11m13L"
    betastat <- m1 %>% filter(term=="regionpfc11m13L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==2){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14r14c11mL)
    namevar = "pfc14r14c11mL"
    betastat <- m1 %>% filter(term=="regionpfc14r14c11mL") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==3){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11L)
    namevar = "pfc11L"
    betastat <- m1 %>% filter(term=="regionpfc11L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==4){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc47L)
    namevar = "pfc47L"
    betastat <- m1 %>% filter(term=="regionpfc47L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==5){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14rc2L)
    namevar = "pfc14rc2L"
    betastat <- m1 %>% filter(term=="regionpfc14rc2L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==6){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14m3225L)
    namevar = "pfc14m3225L"
    betastat <- m1 %>% filter(term=="regionpfc14m3225L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==7){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfcfp10L)
    namevar = "pfcfp10L"
    betastat <- m1 %>% filter(term=="regionpfcfp10L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==8){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfcp3224L)
    namevar = "pfcp3224L"
    betastat <- m1 %>% filter(term=="regionpfcp3224L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==9){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfcrdb10L)
    namevar = "pfcrdb10L"
    betastat <- m1 %>% filter(term=="regionpfcrdb10L") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==10){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14r14c11mR)
    namevar = "pfc14r14c11mR"
    betastat <- m1 %>% filter(term=="regionpfc14r14c11mR") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==11){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11aR)
    namevar = "pfc11aR"
    betastat <- m1 %>% filter(term=="regionpfc11aR") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==12){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfcfp10R)
    namevar = "pfcfp10R"
    betastat <- m1 %>% filter(term=="regionpfcfp10R") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==13){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11b47R)
    namevar = "pfc11b47R"
    betastat <- m1 %>% filter(term=="regionpfc11b47R") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==14){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfcrl10R)
    namevar = "pfcrl10R"
    betastat <- m1 %>% filter(term=="regionpfcrl10R") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==15){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14m3225R)
    namevar = "pfc14m3225R"
    betastat <- m1 %>% filter(term=="regionpfc14m3225R") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==16){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfc2432R)
    namevar = "pfc2432R"
    betastat <- m1 %>% filter(term=="regionpfc2432R") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)
  } else if (k==17){
    tempdf <- tempdf %>% dplyr::rename("pfc"=pfcd10R)
    namevar = "pfcd10R"
    betastat <- m1 %>% filter(term=="regionpfcd10R") %>% select(statistic)
    betastat <- betastat %>% dplyr::rename("betastat"=statistic)    
  }
    
  if (nrow(betastat)==0){
    betastat <- NaN
  } 
  
    mb1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                v_max_wi_lag + v_entropy_wi + pfc)^2 + 
                   rt_lag_sc:last_outcome:pfc +
                   rt_vmax_lag_sc:trial_neg_inv_sc:pfc +
                   (rt_vmax_lag_sc+rt_lag_sc|id)+
                   (1|id:run), tempdf)
  mb1 <- tidy(mb1)
  b2b0 <- mb1 %>% filter(term=="trial_neg_inv_sc:rt_vmax_lag_sc:pfc" 
                         | term=="rt_lag_sc:pfc"
                         | term=="pfc"
                         | term=="rt_lag_sc:last_outcomeOmission:pfc"
                         | term=="v_entropy_wi:pfc"
                         | term=="v_max_wi_lag:pfc"
                         | term=="rt_vmax_lag_sc:pfc"
                         | term=="trial_neg_inv_sc:pfc"
                         | term=="last_outcomeOmission:pfc")
  
  b2b0 <- b2b0 %>% mutate(namevar)
  b2b0 <- b2b0 %>% mutate(betavar)
  b2b0 <- b2b0 %>% mutate(betastat)
  b2b_df <- b2b_df %>% add_row(b2b0)
  
  }
  
  pfc_mdf <- grep(pattern='pfc',colnames(mdf))
  for (q in 1:length(pfc_mdf)){
    tempdf <- mdf %>% select(id,run,rt_csv_sc,trial_neg_inv_sc,rt_lag_sc,rt_vmax_lag_sc,
                            last_outcome,v_max_wi_lag,v_entropy_wi,pfc_mdf[q])
    
    if (q==1){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11m13L)
      namevar = "pfc11m13L"
      betastat <- m1 %>% filter(term=="regionpfc11m13L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==2){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14r14c11mL)
      namevar = "pfc14r14c11mL"
      betastat <- m1 %>% filter(term=="regionpfc14r14c11mL") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==3){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11L)
      namevar = "pfc11L"
      betastat <- m1 %>% filter(term=="regionpfc11L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==4){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc47L)
      namevar = "pfc47L"
      betastat <- m1 %>% filter(term=="regionpfc47L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==5){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14rc2L)
      namevar = "pfc14rc2L"
      betastat <- m1 %>% filter(term=="regionpfc14rc2L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==6){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14m3225L)
      namevar = "pfc14m3225L"
      betastat <- m1 %>% filter(term=="regionpfc14m3225L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==7){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfcfp10L)
      namevar = "pfcfp10L"
      betastat <- m1 %>% filter(term=="regionpfcfp10L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==8){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfcp3224L)
      namevar = "pfcp3224L"
      betastat <- m1 %>% filter(term=="regionpfcp3224L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==9){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfcrdb10L)
      namevar = "pfcrdb10L"
      betastat <- m1 %>% filter(term=="regionpfcrdb10L") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==10){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14r14c11mR)
      namevar = "pfc14r14c11mR"
      betastat <- m1 %>% filter(term=="regionpfc14r14c11mR") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==11){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11aR)
      namevar = "pfc11aR"
      betastat <- m1 %>% filter(term=="regionpfc11aR") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==12){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfcfp10R)
      namevar = "pfcfp10R"
      betastat <- m1 %>% filter(term=="regionpfcfp10R") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==13){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc11b47R)
      namevar = "pfc11b47R"
      betastat <- m1 %>% filter(term=="regionpfc11b47R") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==14){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfcrl10R)
      namevar = "pfcrl10R"
      betastat <- m1 %>% filter(term=="regionpfcrl10R") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==15){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc14m3225R)
      namevar = "pfc14m3225R"
      betastat <- m1 %>% filter(term=="regionpfc14m3225R") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==16){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfc2432R)
      namevar = "pfc2432R"
      betastat <- m1 %>% filter(term=="regionpfc2432R") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)
    } else if (q==17){
      tempdf <- tempdf %>% dplyr::rename("pfc"=pfcd10R)
      namevar = "pfcd10R"
      betastat <- m1 %>% filter(term=="regionpfcd10R") %>% select(statistic)
      betastat <- betastat %>% dplyr::rename("betastat"=statistic)    
    }
    
    if (nrow(betastat)==0){
      betastat <- NaN
    } 
    
    mb1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                v_max_wi_lag + v_entropy_wi + pfc)^2 + 
                   rt_lag_sc:last_outcome:pfc +
                   rt_vmax_lag_sc:trial_neg_inv_sc:pfc +
                   (rt_vmax_lag_sc+rt_lag_sc|id)+
                   (1|id:run), tempdf)
    mb1 <- tidy(mb1)
    b2b0 <- mb1 %>% filter(term=="trial_neg_inv_sc:rt_vmax_lag_sc:pfc" 
                           | term=="rt_lag_sc:pfc"
                           | term=="pfc"
                           | term=="rt_lag_sc:last_outcomeOmission:pfc"
                           | term=="v_entropy_wi:pfc"
                           | term=="v_max_wi_lag:pfc"
                           | term=="rt_vmax_lag_sc:pfc"
                           | term=="trial_neg_inv_sc:pfc"
                           | term=="last_outcomeOmission:pfc")
    
    
    b2b0 <- b2b0 %>% mutate(namevar)
    b2b0 <- b2b0 %>% mutate(betavar)
    b2b0 <- b2b0 %>% mutate(betastat)
    b2b_mdf <- b2b_mdf %>% add_row(b2b0)
  
  }
  
}

save(b2b_df,file="Plot_beta_significance_region_df_corrected")
save(b2b_mdf,file="Plot_beta_significance_region_mdf_corrected")
