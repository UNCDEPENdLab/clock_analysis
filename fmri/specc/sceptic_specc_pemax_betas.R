#entropy change betas
library(readr)
library(lme4)
library(tidyverse)
library(psych)
library(emmeans)
library(hexbin)
library(sjPlot)
base_dir <- "/Users/mnh5174/Box/SCEPTIC_fMRI/specc_betas"
design_dir <- file.path(base_dir, "sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed", "pe_max")
beta_dir <- file.path(base_dir, "with_individual_emo") #for every pairwise emo comparison, these 

setwd(base_dir)

behav <- trial_df <- read.csv("specc_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz") %>%
  mutate(
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
    v_entropy_no5=if_else(run_trial <= 5, NA_real_, v_entropy),
    d_auc_sqrt=if_else(d_auc > 0, NA_real_, sqrt(-1*d_auc)), #only compute the sqrt of d_auc for negative (i.e., reasonable) observations
    v_entropy_sqrt=sqrt(v_entropy),
    rew_om=if_else(score_vba > 0, 1, 0)
  ) %>% #for win/loss maps
  group_by(id, run) %>%
  dplyr::mutate(   #compute rt_swing within run and subject
    rt_csv=rt_csv/1000, #rescale for ease in lmer
    rt_lag = dplyr::lag(rt_csv, 1, order_by=trial),
    rt_vmax_lag = dplyr::lag(rt_vmax, 1, order_by=trial),
    rt_vmax_change = abs(rt_vmax - rt_vmax_lag),
    omission = factor(as.numeric(score_csv > 0), levels=c(0,1), labels=c("Omission", "Reward")),
    omission_lag = lag(omission, order_by=trial),
    v_entropy_lag = dplyr::lag(v_entropy, 1, order_by=trial),
    v_entropy_change = v_entropy - v_entropy_lag, #change in entropy
    v_entropy_change_pos = v_entropy_change*(v_entropy_change > 0),
    v_entropy_change_neg = abs(v_entropy_change*(v_entropy_change < 0)),
    vmax_lag = dplyr::lag(v_max, 1, order_by=trial),
    pe_lag = dplyr::lag(pe_max, 1, order_by=trial),
    ppe_lag = if_else(pe_lag > 0, pe_lag, 0),
    npe_lag = if_else(pe_lag < 0, abs(pe_lag), 0),
    rt_swing = abs( c(NA, diff(rt_csv))),
    rt_swing_sqrt=sqrt(rt_swing)) %>%
  ungroup() %>% dplyr::select(-dataset, -rt_next, -score_next, -clock_onset, -isi_onset, -feedback_onset, -iti_onset, -image) %>%
  mutate(id=factor(id), run=factor(run)) %>%
  group_by(id, run) %>%
  mutate(pe_max_wi=pe_max - mean(pe_max, na.rm=TRUE), pe_max_pm=mean(pe_max, na.rm=TRUE),
         pe_max_wi_lag = dplyr::lag(pe_max_wi, 1, order_by=trial),
         pe_max_wi_lag2 = dplyr::lag(pe_max_wi, 2, order_by=trial),
         v_entropy_change_wi=v_entropy_change - mean(v_entropy_change, na.rm=TRUE), v_entropy_change_pm=mean(v_entropy_change, na.rm=TRUE),
         v_entropy_wi=v_entropy - mean(v_entropy, na.rm=TRUE), v_entropy_pm=mean(v_entropy, na.rm=TRUE)) %>%
  mutate(inv_trial = -100/run_trial, rt_lag_z=as.vector(rt_lag)) %>%
  ungroup()


design <- read.table(file.path(design_dir, "pe_max-Intercept-Age_c-BPD_c-BPD_Age_design.txt"), header=TRUE) %>%
    dplyr::select(feat_input_id, id, SPECC_ID, BPD, Age, Female, Age_c, Age_c, BPD_Age) %>%
  mutate(BPD=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")))

metadata <- read.csv(file.path(beta_dir, "pe_max_cluster_metadata.csv")) %>%
  filter(model=="Intercept-Age_c-BPD_c-BPD_Age") %>% droplevels()
subj_betas <- read.csv(file.path(beta_dir, "pe_max_subj_betas.csv.gz")) %>%
  dplyr::select(-Intercept, -Age, -BPD, -Age_c, -BPD_c, -BPD_Age, -I_Age_c, -BPD_IAge) %>%
  filter(model=="Intercept-Age_c-BPD_c-BPD_Age") %>% droplevels()
run_betas <- read.csv(file.path(design_dir, "pe_max_run_betas.csv.gz")) %>% 
  dplyr::select(-Intercept, -Age, -BPD, -Age_c, -BPD_c, -BPD_Age, -I_Age_c, -BPD_IAge) %>%
  filter(model=="Intercept-Age_c-BPD_c-BPD_Age") %>% droplevels()

#overall betas for pe_max
# sanity checks -- all passed
# s2 <- subj_betas %>% inner_join(design)
# s3 <- s2 %>% inner_join(metadata)
# nrow(s2) == nrow(subj_betas)
# nrow(s3) == nrow(subj_betas)

# subj_betas <- subj_betas %>% inner_join(design) %>% inner_join(metadata) %>%
#   mutate(BPD=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
#          Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")))

#split the target emotion for each emo diff into a separate column
subj_betas <- subj_betas %>% separate(l2_contrast, into=c("l2_contrast", "emotion"), sep="-", fill = "right") %>%
  inner_join(design) %>% inner_join(metadata)


#6 clusters -- the 6th appears to be a bit outside the brain in r adPFC
overall_betas <- subj_betas %>% filter(l2_contrast=="overall" & l3_contrast=="BPD_c") %>%
  dplyr::select(-model, -l1_contrast) %>% droplevels() %>% filter(cluster_number <= 5)

for_fa <- overall_betas %>% dplyr::select(id, cluster_number, cope_value) %>%
  spread(key=cluster_number, value = cope_value, sep="_")

lattice::splom(for_fa)

for_fa_wins <- apply(for_fa %>% dplyr::select(-id), 2, function(col) { psych::winsor(col, trim=0.03 )})

lattice::splom(for_fa_wins)

#principal(for_fa, nfactors = 3)
f1 <- fa(for_fa_wins, nfactors = 1, fm="minres")
f2 <- fa(for_fa_wins, nfactors = 2, fm="minres", rotate="oblimin")
f3 <- fa(for_fa_wins, nfactors = 3, fm="minres", rotate = "oblimin")
f4 <- fa(for_fa_wins, nfactors = 4, fm="minres", rotate = "oblimin")
vss(for_fa_wins, rotate="varimax", fm="ml")

#1-factor solution does quite well -- 61%, and best according to MAP

#for some reason we have 80 for PE analysis -- need to review
f1_scores <- factor.scores(for_fa_wins, f1)$scores %>% as.data.frame() %>% setNames(c("c1_SPL")) %>% bind_cols(for_fa %>% dplyr::select(id))

#join with behavior
behav_w_clus <- behav %>% inner_join(f1_scores) %>% inner_join(design) %>%
  mutate(id=factor(id), run=factor(run)) %>%
  group_by(id, run) %>%
  mutate(pe_max_wi=pe_max - mean(pe_max, na.rm=TRUE), pe_max_pm=mean(pe_max, na.rm=TRUE),
         pe_max_wi_lag = dplyr::lag(pe_max_wi, 1, order_by=trial),
         pe_max_wi_lag2 = dplyr::lag(pe_max_wi, 2, order_by=trial),
         v_entropy_change_wi=v_entropy_change - mean(v_entropy_change, na.rm=TRUE), v_entropy_change_pm=mean(v_entropy_change, na.rm=TRUE),
         v_entropy_wi=v_entropy - mean(v_entropy, na.rm=TRUE), v_entropy_pm=mean(v_entropy, na.rm=TRUE)) %>%
  ungroup()

behav %>% anti_join(f1_scores) %>% pull(id) %>% unique()

summary(m1 <- lmer(rt_csv ~ rt_lag + v_entropy + c1_SPL + (1|id/run), behav_w_clus))

#some action here.
summary(m2 <- lmer(rt_csv ~ rt_lag + v_entropy * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m2, type="3")

summary(m3 <- lmer(rt_csv ~ rt_lag + v_entropy_change * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m3, type="3")

summary(m4 <- lmer(rt_csv ~ rt_lag + v_entropy_change_pm + v_entropy_change_wi * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m4, type="3")

summary(m5 <- lmer(rt_csv ~ rt_lag + v_entropy_pm + v_entropy_wi * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m5, type="3")

plot_model(model=m5, type="pred", terms = c("c1_SPL", "v_entropy_wi", "BPD"))


### Actually get into PE effects

summary(m6 <- lmer(rt_csv ~ rt_lag + pe_max_pm + pe_max_wi_lag * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m6, type="3")

#yep, it's at the level of rt_lag interaction
summary(m7 <- lmer(rt_csv ~ pe_max_pm + rt_lag * pe_max_wi_lag * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m7, type="3")

#add rt_vmax_lag -- actually strengthens the PE effect a bit
summary(m8 <- lmer(rt_csv ~ pe_max_pm + rt_lag * pe_max_wi_lag * c1_SPL * BPD + rt_vmax_lag + (1|id/run), behav_w_clus))
car::Anova(m8, type="3")

#add emo and rewFunc -- strengthens effect even further
summary(m9 <- lmer(rt_csv ~ pe_max_pm + rt_lag * pe_max_wi_lag * c1_SPL * BPD + rt_vmax_lag + emotion + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m9, type="3")

#overall, looks like controls scale their RT stickiness in inverse proportion to SPL
#in BPD, it's largely flat -- so the region does not support sticking with a valuable action
plot_model(m9, type="pred", terms=c("c1_SPL", "rt_lag", "BPD")) #"pe_max_wi_lag", 

#can we use swing as a simpler model? Appears not? We get a strong negative effect of PE, so larger PE+ predicts smaller swing (right)
summary(m8 <- lmer(rt_swing ~ pe_max_pm + pe_max_wi_lag * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m8, type="3")

summary(m8 <- lmer(rt_swing ~ pe_max_pm + pe_max_wi_lag2 * c1_SPL * BPD + (1|id/run), behav_w_clus))
car::Anova(m8, type="3")


summary(m8 <- lmer(rt_swing ~ pe_max_pm + pe_max_wi_lag2 * c1_SPL * BPD + rt_vmax_lag + emotion + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m8, type="3")


####

#Sep2019 analyses: EMO diffs for BPD, including specific emotions for each pairwise

# SCRAM > FEAR
scram_gt_fear <- subj_betas %>% filter(l2_contrast=="scram_gt_fear" & l3_contrast=="BPD_c") %>%
  dplyr::select(-model, -l1_contrast) %>% droplevels() %>% mutate(emotion=if_else(is.na(emotion), "scram_gt_fear", emotion)) %>%
  pivot_wider(names_from=emotion, names_prefix="cope_", values_from=cope_value)
  #spread(key="emotion", value="cope_value")



unique(scram_gt_fear$cluster_number)
unique(scram_gt_fear$label)

for_fa <- scram_gt_fear %>% dplyr::select(id, cluster_number, cope_scram_gt_fear) %>%
  pivot_wider(names_from=cluster_number, values_from = cope_scram_gt_fear, names_prefix = "clus_")

hexplom(for_fa)

for_fa_wins <- apply(for_fa %>% dplyr::select(-id), 2, function(col) { psych::winsor(col, trim=0.03 )})

hexplom(for_fa_wins)

#principal(for_fa, nfactors = 3)
f1 <- fa(for_fa_wins, nfactors = 1, fm="minres")
f2 <- fa(for_fa_wins, nfactors = 2, fm="minres", rotate="oblimin")
f3 <- fa(for_fa_wins, nfactors = 3, fm="minres", rotate = "oblimin")
f4 <- fa(for_fa_wins, nfactors = 4, fm="minres", rotate = "oblimin")
vss(for_fa_wins, rotate="varimax", fm="ml")

f1_scores <- factor.scores(for_fa_wins, f1)$scores %>% as.data.frame() %>% setNames(c("c1_vis")) %>% bind_cols(for_fa %>% dplyr::select(id))

#join with behavior
for_fa_wins <- for_fa %>% mutate_at(vars(starts_with("clus")), ~psych::winsor(., trim=0.03))

behav_w_clus <- behav %>% inner_join(f1_scores) %>% inner_join(design) %>% inner_join(for_fa_wins)

#no meaningful 4-way interactions, thank god
mf2pef2only <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + omission_lag  + c1_vis + BPD + emotion)^3 + 
                      scale(rt_lag):omission_lag:c1_vis +
                      (1|id/run), behav_w_clus)
car::Anova(mf2pef2only, type=3)

#clus 1 in large visual cortex
#clus 2 is cerebellum
#clus 3 is parahippocampal

learnable <- behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV"))

mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_2 + BPD + emotion)^4 + 
                      rt_lag_z:omission_lag:clus_2 +
                      (1|id/run), learnable)
car::Anova(mf2pef2only, type=3)

toplot <- scram_gt_fear %>% inner_join(design) %>%
  gather(key="emotion", value="cope", cope_fear, cope_happy, cope_scram)
#ggplot(for_fa_wins %>% inner_join(design), aes(x=emotion, y=clus_2, color=BPD)) + stat_summary()
ggplot(for_fa_wins %>% inner_join(design), aes(x=BPD, y=clus_2)) + stat_summary()

#pattern:
# - these are clusters with + PE sensitivity to scram, - to fear in HC
#     but opposite pattern in BPD: + or neutral to fear, - to scram
ggplot(toplot, aes(x=BPD, y=cope, color=emotion)) + stat_summary(fun.data=mean_cl_boot, position=position_dodge(width=0.5)) +
  facet_wrap(~cluster_number) + geom_hline(yintercept=0)

pdf("Parahippocampal PE scram fear modulation.pdf", width=6, height=5)
ggplot(toplot %>% filter(cluster_number==3 & emotion !="cope_happy") %>% mutate(cope=cope*1000) %>%
         mutate(emotion=dplyr::recode(emotion, "cope_fear"="fear", "cope_scram"="scrambled")), aes(x=BPD, y=cope, color=emotion)) + 
  stat_summary(fun.data=mean_cl_boot, position=position_dodge(width=0.5), size=1.5) +
  geom_hline(yintercept=0) + theme_bw(base_size=24) + ylab("PE Modulation (a.u.)") +
  scale_color_brewer("Emotion", palette="Dark2") +
  xlab("Group")
  
dev.off()


## Try to bring emotion-specific COPES into LMER
emo_cope <- scram_gt_fear %>% inner_join(design) %>% mutate(id=factor(id)) %>%
  gather(key="emotion", value="cope", cope_fear, cope_happy, cope_scram) %>%
  mutate(row_num=1:n()) %>% select(id, BPD, l2_contrast, l3_contrast, cluster_number, emotion, cope) %>%
  pivot_wider(names_from = cluster_number, names_sep = "clus", values_from = cope, names_prefix = "clus_") %>%
  mutate(emotion=sub("cope_", "", emotion))

behav_emo_cope <- behav %>% inner_join(emo_cope)

#clus 1 vis
mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_1 + BPD + emotion)^4 + 
                      rt_lag_z:omission_lag:clus_1 +
                      (1|id/run), behav_emo_cope)
car::Anova(mf2pef2only, type=3)


#clus 2 cerebellum
mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_2 + BPD + emotion)^4 + 
                      rt_lag_z:omission_lag:clus_2 +
                      (1|id/run), behav_emo_cope)
car::Anova(mf2pef2only, type=3)

#clus 3 parahippo

mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_3 + BPD + emotion)^4 + 
                      rt_lag_z:omission_lag:clus_3 +
                      (1|id/run), behav_emo_cope %>% filter(rewFunc %in% c("IEV", "DEV")))
car::Anova(mf2pef2only, type=3)
mm <- lmer_predict(mf2pef2only, fixat0 = "inv_trial", divide="clus_3")

mf2pef2only <- lmer(rt_csv ~ inv_trial + rt_lag_z * omission_lag * clus_3 * BPD * emotion + 
                      #rt_lag_z:omission_lag:clus_3 +
                      (1|id/run), behav_emo_cope %>% filter(rewFunc %in% c("IEV", "DEV")))
car::Anova(mf2pef2only, type=3)
mm <- lmer_predict(mf2pef2only, fixat0 = "inv_trial", divide="clus_3")


mf2pef2only <- lmer(rt_csv ~ inv_trial + omission_lag + rt_vmax_lag + rt_lag_z * clus_3 * BPD + emotion + 
                      #rt_lag_z:omission_lag:clus_3 +
                      (1|id/run), behav_emo_cope %>% filter(rewFunc %in% c("IEV", "DEV")))
car::Anova(mf2pef2only, type=3)
mm <- lmer_predict(mf2pef2only, fixat0 = c("inv_trial", "rt_vmax_lag"), divide="clus_3")


ggplot(mm %>% filter(emotion != "happy"), aes(x=rt_lag_z, y=rt_csv, color=clus_3, linetype=omission_lag)) + geom_line() + facet_grid(emotion~BPD)


#what does clus 3 do in some simple sense?
mf2pef2only <- lmer(rt_csv ~ inv_trial + rt_vmax_lag + rt_lag_z * clus_3 * BPD + omission_lag + rewFunc + emotion +
                      #rt_lag_z:omission_lag:clus_3 +
                      (1|id/run), behav_emo_cope)# %>% filter(rewFunc %in% c("IEV", "DEV")))
summary(mf2pef2only)
car::Anova(mf2pef2only, type=3)
car::vif(mf2pef2only)

emtrends(mf2pef2only, ~omission_lag, var="clus_3", data=behav_emo_cope)

behav_emo_cope$clus_3 <- behav_emo_cope$clus_3/1000
mf2pef2only <- lmer(rt_csv ~ run_trial + rt_vmax_lag + rt_lag_z * clus_3 + omission_lag + rewFunc + emotion +
                      #rt_lag_z:omission_lag:clus_3 +
                      (1|id/run), behav_emo_cope)# %>% filter(rewFunc %in% c("IEV", "DEV")))
summary(mf2pef2only)
car::Anova(mf2pef2only, type=3)

mm <- lmer_predict(mf2pef2only, fixat0 = c("inv_trial", "rt_vmax_lag"), divide="clus_3")
ggplot(mm %>% filter(emotion != "happy"), aes(x=rt_lag_z, y=rt_csv, color=clus_3, linetype=omission_lag)) + 
  geom_line() + facet_grid(rewFunc~emotion)

pdf("parahippo pe effects.pdf", width=7, height=5)
plot_model(mf2pef2only, type="est",terms=c("clus_3", "rt_lag_z:clus_3"), line.size=1.5, fatten=2) +
  ylab("Unstandardized estimates") +
  ggtitle("") + geom_hline(yintercept = 0, color="black") +
  cowplot::theme_cowplot(font_size = 24)
dev.off()
  
  #"rt_lag_z", 

summary(mm <- lmer(rt_swing ~ rt_vmax_lag + clus_3 * emotion + omission_lag + rewFunc + emotion + run_trial + (1|id/run), behav_emo_cope))
car::Anova(mm, type=3)

mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_3 + BPD + emotion)^4 +  #rt_vmax_lag + 
                      rt_lag_z:omission_lag:clus_3 +
                      (1|id/run), behav_emo_cope %>% filter(rewFunc %in% c("IEV", "DEV")))
car::Anova(mf2pef2only, type=3)
mm <- lmer_predict(mf2pef2only, fixat0 = "inv_trial", divide="clus_3")


# mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_3 + BPD + emotion)^4 + 
#                       rt_lag_z:omission_lag:clus_3 +
#                       (1|id/run), behav_emo_cope %>% filter(emotion != "happy")) #simplify scram > fear
# car::Anova(mf2pef2only, type=3)
# 
# mm <- lmer_predict(mf2pef2only, fixat0 = "inv_trial", divide="clus_3")


#looks like > PE modulation in parahippocampal goes with more stickiness in HC, but not in BPD
ggplot(mm %>% filter(emotion != "happy"), aes(x=rt_lag_z, y=rt_csv, color=clus_3, linetype=omission_lag)) + geom_line() + facet_grid(emotion~BPD)

#ggplot(behav_emo_cope, aes(x=rt_lag_z, y=rt_csv, color=clus_3 > 0, linetype=omission_lag)) + stat_smooth(method="gam") + facet_grid(emotion~BPD)
ggplot(behav_emo_cope %>% filter(rewFunc %in% c("IEV", "DEV")), aes(x=rt_lag_z, y=rt_csv, color=clus_3 > 0, linetype=omission_lag)) + stat_smooth(method="loess") + facet_grid(emotion~BPD)




mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_3 + BPD + emotion)^4 + 
                      rt_lag_z:omission_lag:clus_3 +
                      (1|id/run), behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV")))
car::Anova(mf2pef2only, type=3)

mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag  + clus_1 + BPD + emotion)^4 + 
                      rt_lag_z:omission_lag:clus_1 + rewFunc +
                      (1|id/run), behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV")))
car::Anova(mf2pef2only, type=3)




mf2pef2only <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + pe_max_wi_lag  + clus_2 + BPD + emotion)^3 + 
                      #scale(rt_lag):omission_lag:c1_vis +
                      (1|id/run), behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV")))
car::Anova(mf2pef2only, type=3)

#look at pe x BPD x rt_lag
mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + pe_max_wi_lag + BPD + emotion)^3 + 
                      (1|id/run), behav_w_clus)
car::Anova(mf2pef2only, type = 3)

#pe_max x rt_lag x BPD interaction
#  - Normatively, positive PEs increase sticking with a good response
#  - In BPD, this reinforcement-based stickiness is diminished. PEs don't guide sticking to the same extent
#  - This interaction is significant, but weaker, if we switch to reward/omission only.



#rt_vmax_lag
summary(m3 <- lmer(rt_csv ~ rt_lag + pe_max_wi_lag + c1_vis * BPD * emotion + rt_vmax_lag * BPD * c1_vis * emotion + rewFunc + (1|id/run), 
                   data=learnable))
car::Anova(m3, type="3")

simpler <- lmer(rt_csv ~ rt_lag * inv_trial + rt_vmax_lag * BPD * emotion + rewFunc + (1|id/run), data=learnable)
car::Anova(simpler, type="3")

cm <- lmer_predict(simpler, fixat0 = c("rt_lag", "inv_trial"))
ggplot(cm, aes(x=rt_vmax_lag, y=rt_csv, color=BPD)) + geom_line() + facet_grid(rewFunc~emotion)


emtrends(simpler, ~BPD * emotion, var="rt_vmax_lag", data=learnable)


emtrends(m3, ~BPD * emotion, var="rt_vmax_lag", data=learnable)
library(dependlab)




mf2pef2only <- lmer(rt_csv ~ (inv_trial + rt_lag_z + omission_lag + BPD + emotion)^3 + 
                      (1|id/run), behav_w_clus)
car::Anova(mf2pef2only, type=3)
plot_model(mf2pef2only, type="pred", terms=c("rt_lag_z", "omission_lag", "BPD"))
df <- ggpredict(mf2pef2only, type="fe", terms=c("rt_lag_z", "omission_lag", "BPD"))
ggplot(df, aes(x=x, y=predicted, color=group, linetype=facet)) +geom_line()

#the BPD group shows weaker sensitivity to PEs elongating the response (willingness to wait of sorts)

summary(mf2pef2only)
plot_model(mf2pef2only, type="pred", terms=c("rt_lag_z", "pe_max_wi_lag", "BPD"))
#ggpredict(mf2pef2only, type="fe"), terms=c("pe_max_wi_lag", "BPD"))

#look 

summary(m9 <- lmer(rt_csv ~ rt_lag + pe_max_wi_lag + c1_vis * BPD * emotion + rt_vmax_lag * BPD * c1_vis * emotion + rewFunc + (1|id/run), 
                   data=behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV") & emotion != "happy")))
car::Anova(m9, type="3")




### What about clusters for scram versus fear?

#8 regions
fear_betas <- subj_betas %>% filter(l2_contrast=="scram_gt_fear" & l3_contrast=="BPD_c") %>%
  dplyr::select(-model, -l1_contrast) %>% droplevels() #%>% filter(cluster_number <= 5)

for_fa <- fear_betas %>% dplyr::select(id, cluster_number, cope_value) %>%
  spread(key=cluster_number, value = cope_value, sep="_")

hexplom(for_fa %>% select(-id))

for_fa_wins <- apply(for_fa %>% dplyr::select(-id), 2, function(col) { psych::winsor(col, trim=0.03 )})

hexplom(for_fa_wins)

#principal(for_fa, nfactors = 3)
f1 <- fa(for_fa_wins, nfactors = 1, fm="minres")
f2 <- fa(for_fa_wins, nfactors = 2, fm="minres", rotate="oblimin")
f3 <- fa(for_fa_wins, nfactors = 3, fm="minres", rotate = "oblimin")
f4 <- fa(for_fa_wins, nfactors = 4, fm="minres", rotate = "oblimin")
vss(for_fa_wins, rotate="varimax", fm="ml")

#1-factor solution does reasonably well -- 46%, and best according to MAP

#for some reason we have 80 for PE analysis -- need to review
f1_scores <- factor.scores(for_fa_wins, f1)$scores %>% as.data.frame() %>% setNames(c("c1_vis")) %>% bind_cols(for_fa %>% dplyr::select(id))

#join with behavior
behav_w_clus <- behav %>% inner_join(f1_scores) %>% inner_join(design) %>%
  mutate(group=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")),
         id=factor(id), run=factor(run)
  ) %>%
  group_by(id, run) %>%
  mutate(pe_max_wi=pe_max - mean(pe_max, na.rm=TRUE), pe_max_pm=mean(pe_max, na.rm=TRUE),
         pe_max_wi_lag = dplyr::lag(pe_max_wi, 1, order_by=trial),
         pe_max_wi_lag2 = dplyr::lag(pe_max_wi, 2, order_by=trial),
         v_entropy_change_wi=v_entropy_change - mean(v_entropy_change, na.rm=TRUE), v_entropy_change_pm=mean(v_entropy_change, na.rm=TRUE),
         v_entropy_wi=v_entropy - mean(v_entropy, na.rm=TRUE), v_entropy_pm=mean(v_entropy, na.rm=TRUE)) %>%
  ungroup()


### So, we'd expect action in emo diffs


summary(m9 <- lmer(rt_csv ~ rt_lag + pe_max_wi_lag + c1_vis * BPD * emotion + rt_vmax_lag * BPD * c1_vis * emotion + rewFunc + (1|id/run), 
#                   data=behav_w_clus))
data=behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV") & emotion != "happy")))
car::Anova(m9, type="3")

plot_model(m9, type="pred", terms=c("c1_vis", "emotion", "BPD"))


summary(m9 <- lmer(rt_csv ~ rt_lag + pe_max_wi_lag * c1_vis * BPD * emotion + rt_vmax_lag * BPD * c1_vis * emotion + rewFunc + (1|id/run), 
                   #                   data=behav_w_clus))
                   data=behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV") & emotion != "happy")))
car::Anova(m9, type="3")


summary(m9 <- lmer(rt_csv ~ pe_max_pm + rt_lag * pe_max_wi_lag * c1_vis * BPD * emotion + rt_vmax_lag + emotion + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m9, type="3")

summary(m9 <- lmer(rt_swing ~ pe_max_pm + rt_lag * pe_max_wi_lag * c1_vis * BPD * emotion + rt_vmax_lag + emotion + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m9, type="3")

summary(m9 <- lmer(rt_swing ~ c1_vis * BPD * emotion + rt_vmax_lag + emotion + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m9, type="3")

#what about PE effect on choice quality
npe_only <- filter(behav_w_clus, npe_lag > 0 & rewFunc %in% c("DEV", "IEV"))
summary(negmodel <- lmer(ev ~ emotion*group + rewFunc*group*npe_lag + Age_c*group*emotion + vmax_lag + 
                           npe_lag*group*emotion + (1 | SPECC_ID/run), npe_only))
car::Anova(negmodel, type="3")

summary(negmodel <- lmer(ev ~ emotion*group + rewFunc*group*npe_lag*c1_vis + Age_c*group*emotion + vmax_lag + 
                           npe_lag*group*emotion + (1 | SPECC_ID/run), npe_only))

car::Anova(negmodel, type="3")

plot_model(negmodel, type="pred", terms=c("c1_vis", "group", "npe_lag"))

#the C1_vis cluster does not appear to mediate group differences in the PE shift behavior
test_mobj <- mplusObject(
  TITLE="test",
  VARIABLE="BETWEEN = group c1_vis gxc;
  USEVARIABLES = id ev npe_lag c1_vis group gxc;
  CLUSTER = id;
  WITHIN=npe_lag;",
  DEFINE = "CENTER c1_vis (GRANDMEAN);
  group = group - 1; gxc = group * c1_vis;", #npe_x_c1_vis = npe_lag * c1_vis; 
  ANALYSIS = "TYPE = TWOLEVEL RANDOM;
  ESTIMATOR=BAYES; BITERATIONS=(10000);",
  MODEL = "
  %WITHIN%
  ss | ev ON npe_lag;
  %BETWEEN%
  ss ON group c1_vis gxc;
  c1_vis WITH group;
  ev ON group;",
  OUTPUT="STDYX RESIDUAL;",
  rdata = npe_only %>% filter(emotion=="fear")
)

test_fit <- mplusModeler(test_mobj, modelout = "test_ml.inp", run=TRUE, Mplus_command = "/Users/mnh5174/Applications/Mplus/mplus", hashfilename = FALSE)
screenreg(test_fit$results, single.row=TRUE)


#any effect on the PPE side?
#the C1_vis cluster does not appear to mediate group differences in the PE shift behavior
test_mobj_ppe <- mplusObject(
  TITLE="test",
  VARIABLE="BETWEEN = group c1_vis gxc;
  USEVARIABLES = id ev ppe_lag group gxc c1_vis;
  CLUSTER = id;
  WITHIN=ppe_lag;",
  DEFINE = "CENTER c1_vis (GRANDMEAN);
  group = group - 1; gxc = group * c1_vis;", #ppe_x_c1_vis = ppe_lag * c1_vis; 
  ANALYSIS = "TYPE = TWOLEVEL RANDOM;
  ESTIMATOR=BAYES; BITERATIONS=(10000);",
  MODEL = "
  %WITHIN%
  ss | ev ON ppe_lag;
  %BETWEEN%
  ss ON group c1_vis gxc;
  c1_vis WITH group;
  ev ON group;",
  OUTPUT="STDYX RESIDUAL;",
  rdata = ppe_only %>% filter(emotion=="Fear")
)

test_fit_ppe <- mplusModeler(test_mobj_ppe, modelout = "ppe_bayes.inp", run=TRUE, Mplus_command = "/Users/mnh5174/Applications/Mplus/mplus", hashfilename = FALSE)
screenreg(test_fit_ppe$results, single.row=TRUE)








#### Scram happy

happy_betas <- subj_betas %>% filter(l2_contrast=="scram_gt_happy" & l3_contrast=="BPD_c") %>%
  dplyr::select(-model, -l1_contrast) %>% droplevels() %>% filter(cluster_number != 10) #complex item

for_fa <- happy_betas %>% dplyr::select(id, cluster_number, cope_value) %>%
  spread(key=cluster_number, value = cope_value, sep="_")

hexplom(for_fa)

for_fa_wins <- apply(for_fa %>% dplyr::select(-id), 2, function(col) { psych::winsor(col, trim=0.05 )})

hexplom(for_fa_wins)



#principal(for_fa, nfactors = 3)
f1 <- fa(for_fa_wins, nfactors = 1, fm="minres")
f2 <- fa(for_fa_wins, nfactors = 2, fm="minres", rotate="oblimin")
f3 <- fa(for_fa_wins, nfactors = 3, fm="minres", rotate = "oblimin")
f4 <- fa(for_fa_wins, nfactors = 4, fm="minres", rotate = "oblimin")
vss(for_fa_wins, rotate="varimax", fm="ml")

#accept two factor solution: most everything on F1.
#on F2, we have lateral IPS
f2


f2_scores <- factor.scores(for_fa_wins, f2)$scores %>% as.data.frame() %>% setNames(c("c1_lfp", "c2_ips")) %>% bind_cols(for_fa %>% dplyr::select(id))

#join with behavior
behav_w_clus <- behav %>% inner_join(f2_scores) %>% inner_join(design) %>%
  mutate(group=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")),
         id=factor(id), run=factor(run)
  ) %>%
  group_by(id, run) %>%
  mutate(pe_max_wi=pe_max - mean(pe_max, na.rm=TRUE), pe_max_pm=mean(pe_max, na.rm=TRUE),
         pe_max_wi_lag = dplyr::lag(pe_max_wi, 1, order_by=trial),
         pe_max_wi_lag2 = dplyr::lag(pe_max_wi, 2, order_by=trial),
         v_entropy_change_wi=v_entropy_change - mean(v_entropy_change, na.rm=TRUE), v_entropy_change_pm=mean(v_entropy_change, na.rm=TRUE),
         v_entropy_wi=v_entropy - mean(v_entropy, na.rm=TRUE), v_entropy_pm=mean(v_entropy, na.rm=TRUE)) %>%
  ungroup()


### So, we'd expect action in emo diffs
summary(m9 <- lmer(rt_csv ~ pe_max_pm + rt_lag * pe_max_wi_lag * c1_lfp * BPD * emotion + rt_vmax_lag + emotion + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m9, type="3")


summary(m9 <- lmer(rt_csv ~ rt_lag + pe_max_wi_lag + c1_lfp * BPD * emotion + rt_vmax_lag * BPD * c1_lfp * emotion + rewFunc + (1|id/run), 
                   data=behav_w_clus))
                    #data=behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV") & emotion != "fear")))
car::Anova(m9, type="3")

plot_model(m9, type="pred", terms=c("rt_vmax_lag", "c1_lfp", "BPD"))
plot_model(m9, type="pred", terms=c("emotion", "c1_lfp", "BPD"))


dat <- ggpredict(m9, type="fe", terms=c("rt_vmax_lag", "c1_lfp", "BPD"))
dat$facet <- factor(dat$facet, levels=c(0,1), labels=c("HC", "BPD"))
pdf("PE scram happy bpd betas.pdf", width=8, height=4)
ggplot(dat, aes(x = x/10, y = predicted, colour = group)) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, size=1.4) +
  facet_wrap(~facet) +
  labs(
    colour = "Lateral FP\nbeta",
    y = "Predicted response time (secs)",
    x = "Location of value max on t-1 (secs)",
    title = "Value-guided choices by group and L-FP PE activity"
  ) + scale_color_brewer(palette="Dark2", breaks=c("-1", "0", "1"), labels=c("-1 SD", "M", "+ 1 SD")) +
  theme_bw(base_size=16)
dev.off()




summary(m9 <- lmer(rt_csv ~ rt_lag + pe_max_wi_lag * c1_lfp * BPD + emotion + rt_vmax_lag + rewFunc + (1|id/run), 
#                   data=behav_w_clus))
data=behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV") & emotion != "fear" & pe_max_wi_lag > 0)))
car::Anova(m9, type="3")
plot_model(m9, type="pred", terms=c("pe_max_wi_lag", "c1_lfp", "BPD"))

summary(m9 <- lmer(rt_csv ~ rt_lag + pe_max_wi_lag * rt_vmax * c1_lfp * BPD + emotion + rt_vmax_lag + rewFunc + (1|id/run), 
                   #                   data=behav_w_clus))
                   data=behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV") & emotion != "fear" & pe_max_wi_lag > 0)))


summary(m9 <- lmer(ev ~ rt_lag + pe_max_wi_lag * c1_lfp * BPD + emotion + rt_vmax_lag + rewFunc + (1|id/run), 
                   #                   data=behav_w_clus))
                   data=behav_w_clus %>% filter(rewFunc %in% c("IEV", "DEV") & emotion != "fear")))
car::Anova(m9, type="3")


######



#### ENTROPY LEFTOVERS
summary(m4 <- lmer(rt_csv ~ rt_lag + v_entropy_change_pm + v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m4, type="3") #yep

#add omission lag (doesn't knock out effect)
summary(m4 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + v_entropy_change_pm + 
                     v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m4, type="3") #yep

# bobyqa complains about tolerance, though it's super close to the .002 mark
# behav_w_clus_z <- behav_w_clus %>% mutate_if(is.numeric, list(function(x) { as.vector(scale(x)) }) )
# summary(m4 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + v_entropy_change_pm + v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run), behav_w_clus_z))

emtrends(m4, specs="emotion", var="v_entropy_change_wi") #so, in general, in fear condition, entropy increases go with slow down, and decreases with speed-up.
#but in happy and scram, changes don't affect RTs.

#looks like a cross-over interaction. High values in controls
df <- summary(emtrends(m4, specs= ~ c3_fp*BPD, var="v_entropy_change_wi",
              at=list(c3_fp=quantile(f3_scores$c3_fp, seq(.2, .8, .2)))))

ggplot(df, aes(x=c3_fp, y=v_entropy_change_wi.trend, ymin=v_entropy_change_wi.trend-SE, ymax=v_entropy_change_wi.trend+SE, color=BPD)) + geom_pointrange() + geom_line() +
  geom_rug(data=f3_scores %>% merge(design), aes(x=c3_fp, y=NULL, ymin=NULL, ymax=NULL, color=factor(BPD)))

ggplot(behav_w_clus, aes(y=v_entropy_change_wi, x=c3_fp, color=BPD)) + stat_smooth(method="loess")

clusters <- f3_scores %>% inner_join(design)
#indeed, all factor scores are higher in BPD
t.test(c1_putamen ~ BPD, clusters)
tapply(clusters$c1_putamen, clusters$BPD, mean)
t.test(c2_parahippo ~ BPD, clusters)
tapply(clusters$c2_parahippo, clusters$BPD, mean)
t.test(c3_fp ~ BPD, clusters)
tapply(clusters$c3_fp, clusters$BPD, mean)

#BPD group shows far higher betas in this region, consistent with selection for cluster
ggplot(clusters, aes(x=factor(BPD), y=c3_fp)) + geom_boxplot()


#####
#does this model apply to other clusters?
#parahippo seems to do nothing here
summary(m5 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + v_entropy_change_pm + 
                     v_entropy_change_wi * c2_parahippo * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m5, type="3")


#putamen does nothing for RTs per se-- actually, seems to have some control over stickiness (mirror of rt_swing?)
summary(m6 <- lmer(rt_csv ~ rt_lag*c1_putamen + omission_lag*rewFunc + v_entropy_change_pm + rt_vmax_lag + #rt_lag*c1_putamen*BPD
                     v_entropy_change_wi * c1_putamen * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m6, type="3")

plot_model(m6, type = "pred", terms=c("rt_lag", "c1_putamen")) #, "BPD"))

#deeper dive into putamen
#putamen is controlling magnitude of RT swing -- higher c1 in BPD drives greater RT swings in response to entropy changes

summary(m6 <- lmer(rt_swing ~ omission_lag*rewFunc + v_entropy_change_pm + scale(rt_vmax_lag) +
                     v_entropy_change_wi * c1_putamen * BPD + emotion * v_entropy_change_wi * c1_putamen + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m6, type="3")

# emmip(m6, ~ v_entropy_change_wi * c1_putamen | BPD)
# emtrends(m6, var="v_entropy_change_wi", ~ v_entropy_change_wi * c1_putamen | BPD)

library(sjPlot)

#handy!
plot_model(m6, type = "pred", terms=c("v_entropy_change_wi", "c1_putamen", "BPD"))

#so, in general, in fear condition, entropy increases go with slow down, and decreases with speed-up.
#but in happy and scram, changes don't affect RTs.





##################

run_betas <- run_betas %>% inner_join(design) %>% inner_join(metadata) %>%
  mutate(BPD=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")))

#run betas
overall_betas <- run_betas %>% filter(l2_contrast=="overall" & l3_contrast=="BPD_c") %>%
  dplyr::select(-model, -l1_contrast) %>% droplevels() %>% filter(cluster_number <= 12)

for_fa <- overall_betas %>% dplyr::select(id, run_num, cluster_number, cope_value) %>%
  spread(key=cluster_number, value = cope_value, sep="_") %>% dplyr::select(-id, -run_num)

lattice::splom(for_fa)

vss(for_fa, rotate="varimax", fm="ml")
f2 <- fa(for_fa, nfactors = 2, fm="ml", rotate="oblimin")
f3 <- fa(for_fa, nfactors = 3, fm="ml", rotate = "oblimin")
f4 <- fa(for_fa, nfactors = 4, fm="ml", rotate = "oblimin")

Sx <- cov(for_fa)
D2 <- mahalanobis(for_fa, colMeans(for_fa), Sx)
D2[which(D2 > quantile(D2, .75) + 1.5*IQR(D2))]

for_fa_wins <- apply(for_fa, 2, function(col) { psych::winsor(col, trim=0.02 )})

lattice::splom(for_fa_wins)

vss(for_fa_wins, rotate="varimax", fm="ml")
f2 <- fa(for_fa_wins, nfactors = 2, fm="ml", rotate="oblimin")
f3 <- fa(for_fa_wins, nfactors = 3, fm="ml", rotate = "oblimin")

