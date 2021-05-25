# 2021-03-03 AndyP
# This function plots the results from Plot_beta_significance.R

setwd("~/vmPFC")

library(ggpubr)
library(tidyverse)

load("Plot_beta_significance_region_df_corrected")
load("Plot_beta_significance_region_mdf_corrected")

df <- b2b_df %>% filter((term=="rt_lag_sc:pfc" | term=="rt_vmax_lag_sc:pfc") & (betavar=="rt_vmax_change" | betavar=="v_entropy"))
mdf <- b2b_mdf %>% filter((term=="rt_lag_sc:pfc" | term=="rt_vmax_lag_sc:pfc") & (betavar=="rt_vmax_change" | betavar=="v_entropy"))
#df$signal <-str_split_fixed(df$betavar, "-", 2)[,1]
#df$align <-str_split_fixed(df$betavar, "-", 2)[,2]
df$align <- df$betavar
mdf$align <- mdf$betavar
df$dataset <- "fMRI"
mdf$dataset <- "MEG"

qdf <- df %>% add_row(mdf)

pdf("pfc_beta_behavior_effects4.pdf", height = 5, width = 8)
ggplot(qdf, aes(betastat, statistic, color=as.factor(network),shape=dataset)) + geom_point() + geom_text(size=2,aes(label=namevar), 
                                  position=position_jitter(width=0,height=1)) +  facet_grid(align~term)
#ggplot(qdf, aes(betastat, statistic, color=as.factor(dataset))) + geom_point() + facet_grid(align~term) + stat_cor(aes(betastat,statistic))
dev.off()

tempdf <- df %>% select(statistic)
tempdf <- tempdf %>% rename("dfstat"=statistic)
temp2 <- mdf %>% select(statistic)
tempdf$mdfstat <- temp2$statistic
tempdf$betavar <- df$betavar

pdf("test.pdf",height=5,width=8)
ggplot(tempdf, aes(dfstat,mdfstat, color=betavar)) + geom_point()
dev.off()

PFC_gradients <- read_csv('PFC_region_gradients.csv')

Schaefer_label <- NULL
MLgrad <- NULL
APgrad <- NULL
ISgrad <- NULL
LR <- NULL
symmetry_group <- NULL
network <- NULL
pe_fa <- NULL
v_chosen_fa <- NULL
rt_vmax_change_fa <- NULL
v_entropy_fa <- NULL
nD <- nrow(df)
nR <- nrow(PFC_gradients)
for (i in 1:nD){
  good = FALSE
  for (j in 1:nR){
    if (df$namevar[i] == PFC_gradients$region[j]) {
      Schaefer_label[i] <- PFC_gradients$region[j]
      MLgrad[i] <- PFC_gradients$`M/L`[j]
      APgrad[i] <- PFC_gradients$`A/P`[j]
      ISgrad[i] <- PFC_gradients$`I/S`[j]
      LR[i] <- PFC_gradients$`L/R`[j]
      symmetry_group[i] <- PFC_gradients$symmetry_group[j]
      network[i] <- PFC_gradients$network[j]
      pe_fa[i] <- PFC_gradients$pe_max_fa[j]
      v_chosen_fa[i] <- PFC_gradients$v_chosen_fa[j]
      rt_vmax_change_fa[i] <- PFC_gradients$rt_vmax_change_fa[j]
      v_entropy_fa[i] <- PFC_gradients$v_entropy_fa[j]
      good = TRUE
      break
    }
  }
  if (!any(PFC_gradients$region == df$namevar[i])){
    print(i)
  }
} 

df <- df %>% mutate(MLgrad = MLgrad, APgrad = APgrad, ISgrad = ISgrad, LR = LR, symmetry_group=symmetry_group, network = network,
                    pe_fa = pe_fa, v_chosen_fa = v_chosen_fa, rt_vmax_change_fa = rt_vmax_change_fa, v_entropy_fa = v_entropy_fa)


pdf("test5.pdf", height = 5, width = 8)
ggplot(df, aes(betastat, statistic)) + geom_point() + facet_grid(align~term)+ stat_cor(aes(betastat,statistic))
dev.off()

# factor analysis
#setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/')
#setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-v_chosen-preconvolve_fse_groupfixed/v_chosen/')
#setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-rt_vmax_change-preconvolve_fse_groupfixed/rt_vmax_change/')
setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/')
#Hbetas <- read_csv("pe_max_Schaeffer_2018_vmPFC_mask_betas.csv.gz")
#Hbetas <- read_csv("v_chosen_Schaeffer_2018_vmPFC_mask_betas.csv.gz")
#Hbetas <- read_csv("rt_vmax_change_Schaeffer_2018_vmPFC_mask_betas.csv.gz")
Hbetas <- read_csv("v_entropy_Schaeffer_2018_vmPFC_mask_betas.csv.gz")


h <- as_tibble(Hbetas %>% filter(l2_contrast=="overall") %>% filter(atlas_value %in% c(55,56,65,66,67,84,86,88,89,159,160,161,170,171,191,192,194)))

h <- h %>% mutate(Schaefer_label = case_when(
  atlas_value  == 159 ~ "pfc14r14c11mR",
  atlas_value  == 191 ~ "pfc14m3225R",
  atlas_value  == 160 ~ "pfc11aR",
  atlas_value  == 170 ~ "pfc11b47R",
  atlas_value  == 161 ~ "pfcfp10R",
  atlas_value  == 192 ~ "pfc2432R",
  atlas_value  == 194 ~ "pfcd10R",
  atlas_value  == 171 ~ "pfcl10R",
  atlas_value  == 55 ~ "pfc11m13L",
  atlas_value  == 56 ~ "pfc14r14c11mL",  # limbic
  atlas_value  == 84 ~ "pfc14m3225L",
  atlas_value  == 86 ~ "pfcfp10L",
  atlas_value  == 88 ~ "pfc2432L",
  atlas_value  == 89 ~ "pfcrdb10L",
  atlas_value  == 55 ~ "pfc11m13L",   # limbic
  atlas_value  == 65 ~ "pfc11L",
  atlas_value  == 66 ~ "pfc47L",
  atlas_value  == 67 ~ "pfc14rc2L"
  #atlas_value  == 15 ~ region0[15],
  #atlas_value  == 16 ~ region0[16],
), side = substr(Schaefer_label, nchar(Schaefer_label), nchar(Schaefer_label)),
region_S = substr(Schaefer_label, 1, nchar(Schaefer_label)-1)
)
h_w <- tidyr::pivot_wider(h, names_from=Schaefer_label, values_from=beta, id_cols = id)
just_rois <- h_w[,2:ncol(h_w)]
clust_cor <- cor(just_rois,method = 'pearson',use="complete.obs")
mh <- nfactors(clust_cor, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = TRUE)
h.fa = psych::fa(just_rois, nfactors=3, rotate = "oblimin", fm = "pa") # two factors are vmPFC and central OFC
F <- unclass(h.fa$loadings)
setwd("~/vmPFC")
write.csv(F,"entropy.csv")
