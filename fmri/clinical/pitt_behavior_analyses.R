# brain-to-behavior analyses with anterior and posterior hippocampal cluster betas
# first run explore_beta_cluster_import_pca.R if not run once already

library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(ggpubr)
library(lmerTest)
library(stargazer)
library(car)
library(sjstats)
library(sjPlot)
library(emmeans)
library(cowplot)
library(compareGroups)
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')
# library(stringi)

# source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/clinical/')

explore = T
bsocial = F
plots = F
### load data
 if (explore) {
   load('~/Box/skinner/data/MRI/clock_explore/vba_explore_fixedpara_out.rdata')
   load('~/Box/skinner/data/MRI/clock_explore/explore_subj_df.rdata')
   load('~/Box/skinner/data/MRI/clock_explore/explore_betas.rdata')
   load('~/Box/skinner/data/MRI/clock_explore/explore_allbetas.rdata')
 }
if (bsocial) {
load('~/Box/skinner/data/MRI/clock_bsocial/vba_bsocial_mfx_out.rdata')
load('~/Box/skinner/data/MRI/clock_bsocial/bsocial_subj_df.rdata')
load('~/Box/skinner/data/MRI/clock_bsocial/bsocial_betas.rdata')
load('~/Box/skinner/data/MRI/clock_bsocial/bsocial_allbetas.rdata')

}
# recode run and run_trial
if (bsocial) {
trial_df <- trial_df %>% select(-run, -run_trial) %>%  group_by(id, rewFunc, emotion) %>% mutate(
  run_trial = 1:n(),
  run = floor((trial-1)/n()) + 1
) %>% ungroup()
}

if (explore) {
  trial_df <- trial_df %>% select(-run, -run_trial) %>% mutate(
  run = floor((trial-1)/40) + 1
  )
trial_df <- trial_df %>% group_by(id, run) %>%
  mutate(run_trial = 1:n())
  }


# diagnostics for an incorrect contingency
# ins <- trial_df %>% filter(trial_df$run_trial>50)

# get lags
trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       rt_vmax_lag2 = lag(rt_vmax_lag), # take second lag to decorrelate from rt_lag
                                       v_chosen_lag = lag(v_chosen),
                                       v_max_wi = scale(v_max),
                                       v_max_wi_lag = lag(v_max_wi),
                                       v_entropy_wi = scale(v_entropy),
                                       v_max_b = mean(na.omit(v_max)),
                                       v_entropy_b = mean(na.omit(v_entropy)),
                                       rt_change = rt_csv - rt_lag,
                                       pe_max_lag = lag(pe_max),
                                       abs_pe_max_lag = abs(pe_max_lag),
                                       rt_vmax_change = rt_vmax - rt_vmax_lag,
                                       trial_neg_inv_sc = scale(-1/run_trial),
                                       task_trial_neg_inv_sc = scale(-1/trial),
                                       v_chosen_change = v_chosen - lag(v_chosen)) %>% ungroup() %>%
  mutate(rt_lag_sc = scale(rt_lag),
         rt_csv_sc = scale(rt_csv),
         rt_vmax_lag_sc = scale(rt_vmax_lag),
         rt_vmax_lag2_sc = scale(rt_vmax_lag2),
         rewFunc = fct_relevel(rewFunc, "IEV"))

subject_df$id <- as.character(subject_df$redcapid)
# subject_df$id <- as.character(subject_df$id)
# View(subject_df)
# TEMP: get rid of NA duplicate rows in subject characteristics
subject_df <- subject_df %>%
  arrange(rowSums(is.na(.))) %>%        # sort rows by number of NAs
  distinct(id, .keep_all = TRUE)

if (bsocial){
subject_df <- subject_df %>% mutate(   # recode group
  groupLeth = case_when(
    GroupATT == 'HC' ~ 'HC',
    GroupATT == 'NON' ~ 'BPD_NON',
    GroupATT == 'IDE' ~ 'BPD_NON',
    GroupATT == 'ATT' & Lethality == 'll' ~ 'BPD_LL',
    GroupATT == 'ATT' & Lethality == 'hl' ~ 'BPD_HL'

  )
)
df <- inner_join(trial_df, subject_df) %>% filter(id!=219757 & !is.na(groupLeth) & GroupATT !='89')
df <- inner_join(df, betas)
}
# merge
if (explore) {
  table(subject_df$GroupNEW)
  subject_df <- subject_df %>% mutate(   # recode group
    Group3 = case_when(
      Group == 'HC' ~ 'HC',
      Group == 'ATT' ~ 'ATT',
      Group == 'IDE' ~ 'DNA',
      Group == 'DNA' ~ 'DNA',
      Group == 'DEP' ~ 'DNA'
    ))
  df <- inner_join(trial_df, subject_df) #%>% filter(GroupATT !='88')
  # betas$id <- as.character(betas$ID)
  df <- inner_join(df, betas)
  # lost 2 from subject_df and 4 from trial_df
}
# inspect total earnings by group
sum_df <- trial_df %>% group_by(id) %>% dplyr::summarize(total_earnings = sum(score_csv)) %>% arrange(total_earnings)
if (bsocial) {
sdf <- inner_join(sum_df, subject_df) %>% filter(id!=219757 & !is.na(GroupATT) & GroupATT !='89')
ggplot(sdf, aes(groupLeth, total_earnings)) + geom_boxplot()
summary(lm(total_earnings ~ groupLeth + wtar_raw, sdf))
anova(lm(total_earnings ~ groupLeth, sdf))
sdf <- inner_join(sdf, betas)
asdf <- inner_join(sdf, allbetas)
# no group differences in total earnings
print(pe1 <- createTable(compareGroups(groupLeth ~ pe_max_z_cluster_1  +  pe_max_z_cluster_2 +
                                         pe_max_z_cluster_3   +  pe_max_z_cluster_4   +
                                         pe_max_z_cluster_5 +    pe_max_z_cluster_6 +
                                         pe_max_z_cluster_7   +  pe_max_z_cluster_8   +
                                         pe_max_z_cluster_9 +  pe_max_z_cluster_10  +
                                         pe_max_z_cluster_11, asdf)))
export2html(pe1, "bsocial_clock_pe_clusters_by_group.html")
print(h1 <- createTable(compareGroups(groupLeth ~ v_entropy_z_cluster_1 +
                                        v_entropy_z_cluster_1 + v_entropy_z_cluster_2 + v_entropy_z_cluster_3 +
                                        v_entropy_z_cluster_4 + v_entropy_z_cluster_5 + v_entropy_z_cluster_6 +
                                        v_entropy_z_cluster_7 + v_entropy_z_cluster_8 + v_entropy_z_cluster_9 +
                                        v_entropy_z_cluster_9_overlap +
                                        v_entropy_z_cluster_10 + v_entropy_z_cluster_11, asdf)))
# differences in clusters 8 (R inf. parietal supramarginal gyrus) and 10 (left fusiform/parahippocampal g.)
export2html(h1, "bsocial_clock_h_clusters_by_group.html")

}
if (explore) {
  sdf <- inner_join(sum_df, subject_df) #%>% filter(GroupATT !='88')
  ggplot(sdf, aes(Group, total_earnings)) + geom_boxplot()
  summary(lm(total_earnings ~ GroupNEW + wtar_s + exit_total, sdf))
  sdf <- inner_join(sdf, betas)
  asdf <- inner_join(sdf, allbetas)
  # no group differences in total earnings
  print(pe1 <- createTable(compareGroups(GroupNEW ~ pe_max_cluster_1_3mm +
                                           pe_max_cluster_10_3mm + pe_max_cluster_11_3mm +
                                           pe_max_cluster_1_3mm + pe_max_cluster_2_3mm + pe_max_cluster_3_3mm +
                                           pe_max_cluster_4_3mm + pe_max_cluster_5_3mm + pe_max_cluster_6_3mm +
                                           pe_max_cluster_7_3mm + pe_max_cluster_8_3mm + pe_max_cluster_9_3mm, asdf)))
  export2html(pe1, "explore_clock_pe_clusters_by_group.html")
  print(h1 <- createTable(compareGroups(GroupNEW ~ v_entropy_cluster_1_3mm +
                                           v_entropy_cluster_10_3mm + v_entropy_cluster_11_3mm +
                                           v_entropy_cluster_1_3mm + v_entropy_cluster_2_3mm + v_entropy_cluster_3_3mm +
                                           v_entropy_cluster_4_3mm + v_entropy_cluster_5_3mm + v_entropy_cluster_6_3mm +
                                           v_entropy_cluster_7_3mm + v_entropy_cluster_8_3mm + v_entropy_cluster_9_3mm +
                                          v_entropy_cluster_9_3mm_overlap, asdf)))
  # differences in clusters 8 (R inf. parietal supramarginal gyrus) and 10 (left fusiform/parahippocampal g.)
  export2html(h1, "explore_clock_h_clusters_by_group.html")

  # plot betas by group
  pdf("explore_ph_betas_by_group.pdf", height = 4, width = 4)
  ggplot(sdf, aes(GroupNEW, PH_pe, color = GroupNEW)) + geom_boxplot() + geom_jitter() +
    ylab("Post. hipp. RPE signal") + xlab("") + theme(legend.position = "none")
  dev.off()
  pdf("explore_ah_betas_by_group.pdf", height = 4, width = 4)
  ggplot(sdf, aes(GroupNEW, AH_h_neg_o, color = GroupNEW)) + geom_boxplot() + geom_jitter()+
    ylab("Ant. hipp. global max signal") + xlab("") + theme(legend.position = "none")
  dev.off()
  }

if (explore) {
chars <- sdf %>% select(c(PH_pe, AH_h_neg,AH_h_neg_o, age, contains("total"), contains("raw")))
# careful with the missing data
clust_cor <- corr.test(chars,method = 'pearson', use = "complete.obs")
pdf("explore_hipp_chars_corr.pdf", width=24, height=24)
corrplot(clust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = clust_cor$p, sig.level=0.05, insig = "blank")
dev.off()
}
if (bsocial) {
  # only BPD
  # chars <- sdf  %>% filter(groupLeth!="HC") %>% select(c(PH_pe, AH_h_neg, age, contains("total"), contains("raw"),
                                                         # contains("ipip"), contains("score"), contains("ipde")))
  chars <- sdf %>% select(c(PH_pe, AH_h_neg, age, contains("total"), contains("raw"), contains("ipip"), contains("score")))

    # temporary: delete implausible IPIP_NEO values
  outlier <- function(x) {
    x[x < quantile(x,0.25, na.rm = T) - 6 * IQR(x, na.rm = T) | x > quantile(x,0.75, na.rm = T) + 6 * IQR(x, na.rm = T)] <- NA
    x
  }
  # chars <- chars %>% mutate_if(is.numeric, outlier)
  clust_cor <- corr.test(chars,method = 'pearson', use = "complete.obs")
  pdf("bs_hipp_chars_corr.pdf", width=24, height=24)
  corrplot(clust_cor$r, cl.lim=c(-1,1),
           method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
           order = "hclust", diag = FALSE,
           addCoef.col="black", addCoefasPercent = FALSE,
           p.mat = clust_cor$p, sig.level=0.05, insig = "blank")
  dev.off()
}
# ggplot(sdf %>% filter(sdf$GroupNEW!='HC'), aes(age, PH_pe, color = GroupNEW)) + geom_jitter() + geom_smooth(method = 'gam')
# preliminary group characteristics table
if (bsocial) {
print(c1 <- createTable(compareGroups(groupLeth ~ age + female + edu + wtar_s + exit_total + ipde_cm + spsi_imp_sub +
                                        ipipneo120_o + ipipneo120_c + ipipneo120_e + ipipneo120_a + ipipneo120_n + bis_total, sdf)))
export2html(c1, "bsocial_clock_group_characteristics.html")
}
if (explore) {
  print(c1 <- createTable(compareGroups(GroupNEW ~ age + female + edu + wtar_s + exit_total +
                                          neoffi_neuroticism_total + neoffi_extraversion_total + neoffi_conscientiousness_total +
                                          neoffi_agreeableness_total + total_earnings, sdf),show.n = T))
  export2html(c1, "explore_clock_group_characteristics.html")
  print(c2 <- createTable(compareGroups(GroupNEW ~ age + female + edu + wtar_s + exit_total, sdf),show.n = F))
  export2html(c2, "explore_clock_group_characteristics_brief.html")
  
  # need education
}
# check missingness - a lot in Explore!
library(VIM)
library(mice)
df_aggr = VIM::aggr(sdf, col=mice::mdc(1:2), numbers=TRUE, sortVars=TRUE, labels=names(sdf), cex.axis=.7, gap=3, ylab=c("Proportion of missingness","Missingness Pattern"))

# sanity check on modeling
if (plots) {
ggplot(df, aes(run_trial, rt_vmax, lty = rewFunc)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3))
ggplot(df, aes(run_trial, rt_swing_lr, lty = rewFunc)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3))

ggplot(df, aes(run_trial, pe_max, lty = rewFunc)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3))
ggplot(df, aes(run_trial, v_entropy_wi, lty = rewFunc)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,4))


# ggplot(df, aes(run_trial, rt_csv, lty = rewFunc, color = groupLeth)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3))
# ggplot(df, aes(run_trial, rt_swing_lr, color = groupLeth)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3)) +
#   facet_wrap(~rewFunc)
# ggplot(df, aes(run_trial, score_csv, color = groupLeth)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3)) +
#   facet_wrap(~rewFunc)

ggplot(df, aes(run_trial, rt_csv, color = GroupNEW)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3)) +
  facet_wrap(~rewFunc)
pdf("explore_rt_swings_group_geom_smooth.pdf", height = 4, width = 8)
ggplot(df, aes(run_trial, rt_swing_lr, color = GroupNEW)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,4)) +
  facet_wrap(~rewFunc)
dev.off()
ggplot(df, aes(run_trial, score_csv, color = GroupNEW)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,5)) +
  facet_wrap(~rewFunc)

ggplot(df, aes(run_trial, v_entropy, color = GroupNEW)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,4)) +
  facet_wrap(~rewFunc)

}
### comment out individual graphs for speed
# pdf('inspect_rts_ind.pdf', height = 20, width = 20)
# # ggplot(df, aes(run_trial, rt_csv, color = run)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3)) +
# #   facet_wrap(id~rewFunc)
# ggplot(df, aes(run_trial, rt_csv, color = as.factor(run))) + geom_line() +
#   facet_wrap(id~rewFunc)
# dev.off()

pdf('inspect_rt_swings_ind.pdf', height = 20, width = 20)
ggplot(df, aes(run_trial, rt_swing_lr)) + geom_smooth(method = 'gam',  formula = y~splines::ns(x,3)) +
  facet_wrap(id~rewFunc)
dev.off()
###############
# model-free analyses


# -1/trial does not capture the learning curve perfectly
mf1 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rewFunc +  omission_lag) ^2 +
              (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
summary(mf1)
vif(mf1)
# mf2 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rewFunc + omission_lag + groupLeth)^3 +
mf2 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rewFunc + omission_lag + rt_vmax_lag_sc + GroupNEW)^3 +
              (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
summary(mf2)
Anova(mf2, '3')
vif(mf2)
anova(mf1, mf2)

# just 2-way RTlag*grp interaction
em <- as_tibble(emtrends(mf2,var = "rt_lag_sc", specs = c("GroupNEW", "rewFunc")), data = df %>% filter(!is.na(rt_vmax_lag_sc)))
em$RT_swing = -em$rt_lag_sc.trend
group_labels <- c("Suicide \nattempters", "Depressed, \nhigh-ideation", "Depressed, \nlow-ideation", "Controls")
pdf("explore_rt_swings_group_3way.pdf", height = 4, width = 6)
ggplot(em, aes(GroupNEW, RT_swing, color = GroupNEW, lty = rewFunc)) + xlab("") +
  geom_errorbar(aes( ymin = -asymp.LCL, ymax = -asymp.UCL), size = 1) + geom_point(size = 4) + 
  scale_x_discrete(labels = group_labels)  + theme(legend.title = element_blank())
dev.off()

em <- as_tibble(emtrends(mf2,var = "rt_lag_sc", specs = c("GroupNEW", "omission_lag", "rewFunc")), data = df %>% filter(!is.na(rt_vmax_lag_sc)))
em$RT_swing = -em$rt_lag_sc.trend
em$reward <- "Reward"
em$reward[em$omission_lag] <- "Omission"
pdf("explore_rt_swings_group_rewFunc.pdf", height = 4, width = 9)
ggplot(em, aes(GroupNEW, RT_swing, color = GroupNEW, lty = reward)) + xlab("") +
  geom_errorbar(aes( ymin = -asymp.LCL, ymax = -asymp.UCL)) + geom_point() + facet_wrap(~rewFunc)
dev.off()
# control for exit and age
mf2ae <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rewFunc + omission_lag + groupLeth)^3 +
# mf2ae <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rewFunc + omission_lag + GroupNEW)^3 +
                (trial_neg_inv_sc + rt_lag_sc + rewFunc + omission_lag + scale(age))^3 +
                (trial_neg_inv_sc + rt_lag_sc + rewFunc + omission_lag + scale(exit_raw))^3 +
                (1|id/run), df)
anova(mf2ae, '3')
summary(mf2ae)
screen.lmerTest(mf2ae)
ggplot(df, aes(run_trial, rt_swing_lr, color = GroupNEW, lty = rewFunc)) + geom_smooth()
ggplot(df, aes(run_trial, rt_swing, color = GroupNEW, lty = rewFunc)) + geom_smooth()
ggplot(df, aes(run_trial, rt_swing_lr, color = groupLeth)) + geom_smooth(method = "gam") + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, rt_swing, color = groupLeth)) + geom_smooth(method = "gam")  + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, rt_csv, color = groupLeth)) + geom_smooth(method = "gam") + facet_wrap(~rewFunc)


# Model-based analyses -- no need to include 3-way interactions between design variables
# preliminary models for Explore
emb1 <- lmer(rt_csv_sc ~ (rt_lag_sc + rt_vmax_lag_sc + omission_lag) ^2 +
              (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
summary(emb1)

# emb2 <- lmer(rt_csv_sc ~ (rt_lag_sc + rt_vmax_lag_sc + omission_lag + PH_pe) ^2 +
#                (rt_lag_sc + rt_vmax_lag_sc + omission_lag + AH_h_neg) ^2 +
#                rt_lag_sc:omission_lag:PH_pe +
#                rt_lag_sc:omission_lag:PH_pe +
#                (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
# summary(emb2)

embo2 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + PH_pe) ^2 +
               (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + AH_h_neg_o) ^2 +
                rt_lag_sc:omission_lag:PH_pe +
                rt_lag_sc:omission_lag:AH_h_neg_o +
                rt_vmax_lag_sc:trial_neg_inv_sc:AH_h_neg_o +
                rt_vmax_lag_sc:trial_neg_inv_sc:PH_pe +
                (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
summary(embo2)
Anova(embo2, '3')
#
# # odd result in bsocial, is it group?
# embo3 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + PH_pe + groupLeth) ^3 +
#                 (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + AH_h_neg + groupLeth) ^3 +
#                 rt_lag_sc:omission_lag:PH_pe +
#                 # rt_lag_sc:omission_lag:AH_h_neg_o +
#                 rt_vmax_lag_sc:trial_neg_inv_sc:AH_h_neg +
#                 # rt_vmax_lag_sc:trial_neg_inv_sc:PH_pe +
#                 rt_lag_sc:omission_lag:PH_pe:groupLeth +
#                 # rt_lag_sc:omission_lag:AH_h_neg_o:groupLeth +
#                 rt_vmax_lag_sc:AH_h_neg:groupLeth +
#                 # rt_vmax_lag_sc:trial_neg_inv_sc:PH_pe:groupLeth +
#                 (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
embo3 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + PH_pe + GroupNEW) ^3 +
                (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + AH_h_neg_o + GroupNEW) ^3 +
                rt_lag_sc:omission_lag:PH_pe +
                # rt_lag_sc:omission_lag:AH_h_neg_o +
                rt_vmax_lag_sc:trial_neg_inv_sc:AH_h_neg_o +
                # rt_vmax_lag_sc:trial_neg_inv_sc:PH_pe +
                rt_lag_sc:omission_lag:PH_pe:GroupNEW +
                # rt_lag_sc:omission_lag:AH_h_neg_o:GroupNEW +
                rt_vmax_lag_sc:AH_h_neg_o:GroupNEW +
                # rt_vmax_lag_sc:trial_neg_inv_sc:PH_pe:GroupNEW +
                (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))


summary(embo3)
screen.lmerTest(embo3, .01)
Anova(embo3, '3')
vif(embo3)
anova(embo2, embo3)

# bsocial verions
# em <- as_tibble(emtrends(embo3,var = "rt_lag_sc", specs = c("groupLeth", "PH_pe", "omission_lag"), at = list(PH_pe = c(-2,0,2))))
# em$RT_swing = -em$rt_lag_sc.trend
# ggplot(em, aes(PH_pe, RT_swing, color = omission_lag)) + geom_errorbar(aes( ymin = -asymp.LCL, ymax = -asymp.UCL)) + geom_point() + facet_wrap(~groupLeth )
#
# em <- as_tibble(emtrends(embo3,var = "rt_lag_sc", specs = c("groupLeth", "AH_h_neg", "omission_lag"), at = list(AH_h_neg = c(-2,0,2))))
# em$RT_swing = -em$rt_lag_sc.trend
# ggplot(em, aes(AH_h_neg, RT_swing, color = omission_lag)) + geom_errorbar(aes( ymin = -asymp.LCL, ymax = -asymp.UCL)) + geom_point() + facet_wrap(~groupLeth )

# explore versions
em <- as_tibble(emtrends(embo3,var = "rt_lag_sc", specs = c("GroupNEW", "PH_pe", "omission_lag"), at = list(PH_pe = c(-1.21,2.07)), data = df %>% filter(!is.na(rt_vmax_lag_sc))))
em$RT_swing = -em$rt_lag_sc.trend
em$reward <- "Reward"
em$reward[em$omission_lag] <- "Omission"
pdf("explore_PH_RTswings_group.pdf", height = 6, width = 6)
ggplot(em, aes(PH_pe, RT_swing, color = GroupNEW, lty = reward)) + geom_errorbar(aes( ymin = -asymp.LCL, ymax = -asymp.UCL)) + geom_point() +
  facet_wrap(~GroupNEW )
dev.off()
em <- as_tibble(emtrends(embo3,var = "rt_vmax_lag_sc", specs = c("GroupNEW", "AH_h_neg_o"), at = list(AH_h_neg_o = c(-1.68,2.16)), data = df %>% filter(!is.na(rt_vmax_lag_sc))))
em$RTvmax_effect = em$rt_vmax_lag_sc.trend
ggplot(em, aes(AH_h_neg_o, RTvmax_effect, color = GroupNEW)) + geom_errorbar(aes( ymin = asymp.LCL, ymax = asymp.UCL)) + geom_point() + facet_wrap(~GroupNEW )


# just group
embo4 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_entropy_wi  + GroupNEW) ^3 +
                (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
summary(embo4)
screen.lmerTest(embo4)
Anova(embo4, '3')


# compare entropy timecourses across groups

hm1 <- lmer(v_entropy ~ (trial_neg_inv_sc + rewFunc + GroupNEW) ^3 +
              (1|id/run), df)
summary(hm1)
Anova(hm1, '3')
###########
mb1 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag2_sc + omission_lag) ^2 +
              (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
summary(mb1)
anova(mf1, mb1)
vif(mf1)

# mb1a <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_entropy_wi + v_entropy_b) ^2 +
#               (1|id/run), df %>% filter(!is.na(rt_vmax_lag_sc)))
# summary(mb1a)
# anova(mb1, mb1a)

mb2 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag2_sc + omission_lag) ^2 +
              rt_lag_sc*omission_lag*groupLeth  + omission_lag*rt_vmax_lag2_sc*groupLeth +  (1|id/run), df)
summary(mb2) # RT swings: HC=HL < LL=NON
Anova(mb2, '3')
em <- as_tibble(emtrends(mb2, var = 'rt_vmax_lag2_sc', specs = c('omission_lag', 'groupLeth')))
em$reward <- 'Reward'
em$reward[em$omission_lag] <- 'Omission'
ggplot(em, aes(reward, rt_vmax_lag2_sc.trend, color = groupLeth)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = rt_vmax_lag2_sc.trend - SE, ymax = rt_vmax_lag2_sc.trend + SE))

em <- as_tibble(emtrends(mb2, var = 'rt_lag_sc', specs = c('omission_lag', 'groupLeth')))
em$reward <- 'Reward'
em$reward[em$omission_lag] <- 'Omission'
ggplot(em, aes(reward, rt_lag_sc.trend, color = groupLeth)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = rt_lag_sc.trend - SE, ymax = rt_lag_sc.trend + SE))

# control for age and WTAR
mb2aw <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag) ^2 +
                rt_lag_sc*omission_lag*groupLeth  + omission_lag*rt_vmax_lag_sc*groupLeth  +
                rt_lag_sc*omission_lag*scale(age) + omission_lag*rt_vmax_lag_sc*scale(age) +
                rt_lag_sc*omission_lag*scale(wtar_raw)  + omission_lag*rt_vmax_lag_sc*scale(wtar_raw) +
              (1|id/run), df)
summary(mb2aw) # RT swings: HC=HL < LL=NON
screen.lmerTest(mb2aw)

mb2ae <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag) ^2 +
                rt_lag_sc*omission_lag*groupLeth  + omission_lag*rt_vmax_lag_sc*groupLeth  +
                rt_lag_sc*omission_lag*scale(age) + omission_lag*rt_vmax_lag_sc*scale(age) +
                rt_lag_sc*omission_lag*scale(exit_raw)  + omission_lag*rt_vmax_lag_sc*scale(exit_raw) +
                (1|id/run), df)
summary(mb2ae)
screen.lmerTest(mb2ae)
#
save(file = 'bs_models.Rdata', list = ls(all.names = TRUE))
# load('bs_models.Rdata')
