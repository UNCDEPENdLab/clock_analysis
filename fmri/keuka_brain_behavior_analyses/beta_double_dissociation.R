# double dissociation of betas along the long axis
library(wesanderson)
library(modelr)
library(tidyverse)
library(lme4)
library(multcomp)
library(afex)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)
library(viridis)
library(emmeans)
library(RColorBrewer)
library(nlme)

pal = wes_palette("Zissou1", 12, type = "continuous")
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/h_pe_betas')

# make long dataframe of both beta flavors
#selective maintenance

#zstat
h <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_atlas_zstats.csv.gz') %>%
  filter(l2_contrast=="overall") %>% mutate(zstat=-1*zstat) %>% dplyr::rename(beta=zstat)
#range(h$zstat)
#mutate(beta=-1*beta) 
#mutate(beta = -1*scale(beta, center = T))

#model that contains both entropy and pe in L1 FSL model
h_simult <- read_csv('~/Box/SCEPTIC_fMRI/1h_2h_analyses/entropy_pemax_simult/v_entropy_atlas_zstats.csv.gz') %>%
  filter(l2_contrast=="overall") %>% mutate(zstat=-1*zstat) %>% dplyr::rename(beta=zstat)

#raw betas
# h <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_atlas_betas.csv.gz') %>% filter(l2_contrast=="overall") %>%
#    mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
#    group_by(numid) %>% mutate(beta = -1*scale(beta, center = F)) %>% ungroup() #per subject rms transformation

#mutate(beta = -1*scale(beta, center = F)) %>% 


# pdf("entropy_hist_by_bin_selective.pdf", width=12, height=12)
# h$axis_bin <- cut(h$atlas_value, breaks = 12, labels = 1:12)
# lattice::histogram(~zstat | axis_bin, h)
# dev.off()


pe <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_betas.csv.gz") %>%
  filter(l2_contrast=="overall") %>%
  mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(numid) %>% mutate(beta = scale(beta, center = F)) %>% ungroup() #per subject rms transformation

pe <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_zstats.csv.gz') %>%
  filter(l2_contrast=="overall") %>% mutate(zstat=psych::winsor(zstat, trim=0.001)) %>% dplyr::rename(beta=zstat)

#model that also includes entropy
pe_simult <- read_csv('~/Box/SCEPTIC_fMRI/1h_2h_analyses/entropy_pemax_simult/pe_max_atlas_zstats.csv.gz') %>%
  filter(l2_contrast=="overall") %>% mutate(zstat=psych::winsor(zstat, trim=0.001)) %>% dplyr::rename(beta=zstat)


# hist(pe$zstat)
# pe <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_betas.csv.gz") %>%
#   filter(l2_contrast=="overall") %>%
#   mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
#   group_by(numid) %>% mutate(beta = scale(beta, center = F)) %>% ungroup() #per subject rms transformation

  #mutate(beta = scale(beta, center = F))

#full maintenance Fixed LR V

#h_full <- h_full %>% filter(l2_contrast=="overall") %>% mutate(l1_contrast="v_entropy_fixed", beta=-1*beta) %>% rename(ID = id) #beta = -scale(beta, center = F)
#h_full <- read_csv('v_entropy_atlas_zstats.csv.gz')
#h_full <- h_full %>% filter(l2_contrast=="overall") %>%  mutate(l1_contrast="v_entropy_fixed", zstat=-1*zstat)


h_full <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fixed_groupfixed/v_entropy/v_entropy_atlas_zstats.csv.gz') %>%
  filter(l2_contrast=="overall") %>%
  mutate(zstat=-1*zstat, l1_contrast="v_entropy_fixed") %>% dplyr::rename(beta=zstat)

  #mutate(zstat=-1*psych::winsor(zstat, trim=0.001), l1_contrast="v_entropy_fixed") %>% dplyr::rename(beta=zstat)
# h_full <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fixed_groupfixed/v_entropy/v_entropy_atlas_betas.csv.gz') %>%
#  mutate(beta=psych::winsor(beta, trim=0.01)) %>% #1% winsorize before normalization
#  dplyr::rename(ID = id) %>%
#  group_by(numid) %>% mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = F))) %>% ungroup()
#mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = F))) %>%


#mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = F))) %>% rename(ID = id)
#mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = T)))

pe_fixed <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_trial_fixed-preconvolve_fse_groupfixed/pe_trial_fixed/pe_trial_fixed_atlas_betas.csv.gz")
pe_fixed <- pe_fixed %>% filter(l2_contrast=="overall") %>% 
  rename(ID = id) %>%
  mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(numid) %>% mutate(beta = as.vector(scale(beta, center = F))) %>% ungroup()
  
#Revision 1st half/2nd half split betas
h_1h <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/v_entropy_1h_atlas_zstats.csv.gz") %>%
  filter(l2_contrast=="overall") %>% mutate(zstat=-1*zstat) %>% dplyr::rename(beta=zstat)

h_2h <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/v_entropy_2h_atlas_zstats.csv.gz") %>%
  filter(l2_contrast=="overall") %>% mutate(zstat=-1*zstat) %>% dplyr::rename(beta=zstat)

pe_1h <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/pe_1h_atlas_zstats.csv.gz") %>%
  filter(l2_contrast=="overall") %>% dplyr::rename(beta=zstat)

pe_2h <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/pe_2h_atlas_zstats.csv.gz") %>%
  filter(l2_contrast=="overall") %>% dplyr::rename(beta=zstat)

hist(h_1h$beta)
hist(h_2h$beta)
hist(pe_1h$beta)
hist(pe_2h$beta)
xtabs(~atlas_name + l2_contrast, h_1h)
xtabs(~atlas_name + l2_contrast, h_2h)
xtabs(~atlas_name + l2_contrast, pe_1h)
xtabs(~atlas_name + l2_contrast, pe_2h)

df_split <- rbind(h_1h, h_2h, pe_1h, pe_2h)
rm(h_1h, h_2h, pe_1h, pe_2h)

#first-level z stat approach
#pe trial-fixed learning rates .05, .10, .15, .20
# pe_05 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p05_atlas_zstats.csv.gz") %>%
#   filter(l2_contrast=="overall") %>% mutate(zstat=psych::winsor(zstat, trim=0.001)) %>% dplyr::rename(beta=zstat)
# 
# pe_10 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p10_atlas_zstats.csv.gz") %>%
#   filter(l2_contrast=="overall") %>% mutate(zstat=psych::winsor(zstat, trim=0.001)) %>% dplyr::rename(beta=zstat)
# 
# pe_15 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p15_atlas_zstats.csv.gz") %>%
#   filter(l2_contrast=="overall") %>% mutate(zstat=psych::winsor(zstat, trim=0.001)) %>% dplyr::rename(beta=zstat)
# 
# pe_20 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p20_atlas_zstats.csv.gz") %>%
#   filter(l2_contrast=="overall") %>% mutate(zstat=psych::winsor(zstat, trim=0.001)) %>% dplyr::rename(beta=zstat)

#RMS beta approach
pe_05 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p05_atlas_betas.csv.gz") %>%
  filter(l2_contrast=="overall") %>% rename(ID = id) %>% mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(numid) %>% mutate(beta = scale(beta, center = F)) %>% ungroup() #per subject rms transformation

pe_10 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p10_atlas_betas.csv.gz") %>%
  filter(l2_contrast=="overall") %>% rename(ID = id) %>% mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(numid) %>% mutate(beta = scale(beta, center = F)) %>% ungroup() #per subject rms transformation

pe_15 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p15_atlas_betas.csv.gz") %>%
  filter(l2_contrast=="overall") %>% rename(ID = id) %>% mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(numid) %>% mutate(beta = scale(beta, center = F)) %>% ungroup() #per subject rms transformation

pe_20 <- read_csv("~/Box/SCEPTIC_fMRI/1h_2h_analyses/trial_fixed_lr/pe_trial_fixed_p20_atlas_betas.csv.gz") %>%
  filter(l2_contrast=="overall") %>% rename(ID = id) %>% mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(numid) %>% mutate(beta = scale(beta, center = F)) %>% ungroup() #per subject rms transformation

pe_compare_df <- rbind(pe, pe_05, pe_10, pe_15, pe_20)

hdf <- rbind(h, h_full) #%>% dplyr::rename(beta=zstat)
df <- rbind(h, pe)
adf <- rbind(df, pe_fixed)
pedf <- rbind(pe, pe_fixed)

df_simult <- rbind(h_simult, pe_simult)

# there is a heteroscedasticity with greater variance anteriorly
pdf('h_pe_betas_scatter.pdf', height = 6, width = 7)
ggplot(df, aes(atlas_value,beta, color = l1_contrast)) + geom_jitter(alpha = .4)
dev.off()

# pdf('h_pe_betas_line.pdf', height = 6, width = 7)
# ggplot(df, aes(atlas_value,beta, color = l1_contrast)) + geom_line(aes(group = c(numid, l1_contrast)), alpha = .2)
# dev.off()

pdf('h_pe_betas_boxplot.pdf', height = 6, width = 20)
ggplot(df, aes(as.factor(atlas_value),beta, color = l1_contrast)) + geom_boxplot()
dev.off()

pdf('h_pe_betas_smooth.pdf', height = 6, width = 7)
ggplot(df, aes(atlas_value,beta, color = l1_contrast)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,8))
dev.off()

pdf('h_pe_petrialfixed_betas_smooth.pdf', height = 6, width = 7)
ggplot(adf, aes(atlas_value,beta, color = l1_contrast)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,8))
dev.off()
#########
# maps - I think we need z stats to weigh the betas or z stats instead of betas
h1 <- h %>% mutate(low_h_beta = -winsor(beta, trim = .01))
ggplot(h1, aes(y, z, fill = low_h_beta)) + geom_tile() + scale_fill_viridis_c()
pdf('individual_low_h_hippo_maps_xy.pdf', height = 10, width = 10)
ggplot(h1, aes(x, y, fill = low_h_beta)) + geom_tile() + scale_fill_viridis_c() + facet_wrap(~ID)
dev.off()

pe1 <- pe %>% mutate(pe_beta = winsor(beta, trim = .01))
ggplot(pe1, aes(y, z, fill = pe_beta)) + geom_tile() + scale_fill_viridis_c() + facet_wrap(~ID)
pdf('individual_pe_hippo_maps_xy.pdf', height = 10, width = 10)
ggplot(pe1, aes(x, y, fill = pe_beta)) + geom_tile() + scale_fill_viridis_c() + facet_wrap(~ID)
dev.off()

# what about the difference map?
diff <- pe %>% mutate(pe_beta = scale(winsor(beta, .01)),
                      h_beta = -scale(winsor(h$beta, .01)),
                      beta_diff_h_minus_pe = h_beta - pe_beta)
ggplot(diff, aes(x, y, fill = beta_diff_h_minus_pe)) + geom_tile() + scale_fill_viridis_c() + facet_wrap(~ID)
ggplot(diff, aes(z, y, fill = beta_diff_h_minus_pe)) + geom_tile() + scale_fill_viridis_c()
# formal models

# start with linear location
m1 <- lmer(beta ~ atlas_value * l1_contrast + (1|numid), df)
summary(m1)
Anova(m1, '3')

# completely general location
# 12 bins
df$axis_bin <- cut(df$atlas_value, breaks = 12, labels = 1:12)
pedf$axis_bin <- cut(pedf$atlas_value, breaks = 12, labels = 1:12)

m2 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), df)
summary(m2)
Anova(m2, '3')
em <- as_tibble(emmeans(m2, ~ axis_bin | l1_contrast))
em$beta <- em$emmean
# going a little crazy on the details here
pal = wes_palette("Zissou1", 12, type = "discrete")
pdf('h_pe_betas_general_location_12_bins.pdf', height = 3, width = 4.5)
ggplot(em, aes(axis_bin, beta, lty = l1_contrast, color = as.numeric(axis_bin), shape = l1_contrast)) + geom_point(position = position_dodge(width = .5), size = 2.5) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = .5), size = .5) + 
  # scale_color_viridis_d(option = "plasma", guide = 'none') +
  scale_color_gradientn(colors = pal, guide = 'none') +
  # scale_color_brewer(type = 'seq', guide = 'none') +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'), 
        axis.text.x = element_blank())  +
  scale_linetype(labels = c("Prediction error", "Entropy, reversed")) +   scale_shape(labels = c("Prediction error", "Entropy, reversed")) + 
  xlab("Post. <= Long axis position => Ant.")
dev.off()
# 100 bins
df$location <- as.factor(round(df$atlas_value, digits = 2))
pedf$location <- as.factor(round(pedf$atlas_value, digits = 2))

m3 <- lmer(beta ~ location * l1_contrast + (1|numid), df)
summary(m3)
Anova(m3)
em <- as_tibble(emmeans(m3, ~ location | l1_contrast))
em$beta <- em$emmean
ggplot(em, aes(location, beta, color = l1_contrast)) + geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))

# comparisons between SCEPTIC and Q-learning PEs
m4 <- lmer(beta ~ atlas_value * l1_contrast + (1|numid), pedf)
summary(m4)
Anova(m4, '3')

# completely general location
# 12 bins
# pe full versus selective, as reported in paper
m5 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), pedf)
summary(m5)
Anova(m5, '3')

res <- pairs(emmeans(m5, ~  l1_contrast | axis_bin))
summary(glht(m5, linfct=res@linfct), test=adjusted("single-step")) #all pairs at once -- slices 1-3 survive

# 100 bins
df$location <- as.factor(round(df$atlas_value, digits = 2))
m6 <- lmer(beta ~ location * l1_contrast + (1|numid), df)
summary(m6)
Anova(m6)
 em <- as_tibble(emmeans(m3, ~ location | l1_contrast))
em$beta <- em$emmean
ggplot(em, aes(location, beta, color = l1_contrast)) + geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))


#Entropy full versus selective
# completely general location
# 12 bins
hdf$axis_bin <- cut(hdf$atlas_value, breaks = 12, labels = 1:12)
hdf$isfull <- as.numeric(hdf$l1_contrast=="v_entropy_fixed")
hdf$isselective <- as.numeric(hdf$l1_contrast=="v_entropy")

#m5 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), hdf)
m5 <- lmer(beta ~ axis_bin * l1_contrast + (0+isfull|numid) + (0+isselective|numid), hdf)
summary(m5)
Anova(m5, '3')

summary(emmeans(m5, ~ axis_bin | l1_contrast), infer=TRUE)

#does not improve
# m6 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid) + (1|axis_bin), hdf)
# summary(m6)
# 
# anova(m5, m6)

em <- as_tibble(emmeans(m5, ~ axis_bin | l1_contrast))
em$beta <- em$emmean
# going a little crazy on the details here


en_full_sel <- ggplot(em, aes(axis_bin, beta, color = as.numeric(axis_bin), shape = l1_contrast)) + 
  geom_hline(yintercept = 0) +
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey85", show.legend = FALSE) + 
  geom_linerange(aes(ymin = beta-SE, ymax = beta+SE, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .8, color="grey80", show.legend = FALSE) + 
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=as.numeric(axis_bin), linetype=NULL), position = position_dodge(width = .8), size = .4) + 
  #geom_line(aes(x=as.numeric(axis_bin), shape=NULL, linetype=l1_contrast, color = NULL), position = position_dodge(width = .8), size=1) + 
  geom_line(aes(x=as.numeric(axis_bin), shape=NULL, linetype=NULL, color = NULL, group=l1_contrast), position = position_dodge(width = .8), size=1) + 
  geom_point(color="grey10", position = position_dodge(width = .8), size = 5) + 
  geom_point(position = position_dodge(width = .8), size = 4) + 
  
  # scale_color_viridis_d(option = "plasma", guide = 'none') +
  scale_color_gradientn(colors = pal, guide = 'none') +
  # scale_color_brewer(type = 'seq', guide = 'none') +
  theme_gray(base_size=14) + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey65"), 
        panel.grid.minor = element_line(colour = "grey65"), 
        panel.background = element_rect(fill = 'grey60'), 
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2.5 , "lines"),
        legend.text = element_text(lineheight = 1.0),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=6)))  +
  scale_linetype(labels = c("Entropy\nselective", "Entropy\nfull")) + 
  scale_shape(labels = c("Entropy\nselective", "Entropy\nfull")) + 
  #xlab("Post. <= Long axis position => Ant.") +
  xlab("Long axis position") + ylab("Response magnitude (AU)") +
  coord_cartesian(xlim=c(0.8,12.2))

ggsave('h_full_selective_betas_general_location_12_bins_zstat.pdf', en_full_sel, height = 3, width = 5, useDingbats=FALSE)
#ggsave('h_full_selective_betas_general_location_12_bins_persubj_rms.pdf', height = 3, width = 5, useDingbats=FALSE)


## PE trial-wise versus state-wise
pedf$axis_bin <- cut(pedf$atlas_value, breaks = 12, labels = 1:12)
pedf$isfixed <- as.numeric(pedf$l1_contrast=="pe_trial_fixed")
pedf$issceptic <- as.numeric(pedf$l1_contrast=="pe_max")

m6 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), pedf)
summary(m6)
Anova(m6, '3')

emm <- pairs(emmeans(m6, ~l1_contrast | axis_bin)) #, infer=TRUE, adjust="asymptotic")
summary(emm, infer=TRUE, adjust="asymptotic")

summary(glht(m6, linfct=emm@linfct), test=adjusted("single-step"))

m7 <- lmer(beta ~ axis_bin * l1_contrast + (0 + isfixed|numid) + (0 + issceptic | numid), pedf)
summary(m7)
Anova(m7, '3')

anova(m6, m7)

em <- as_tibble(emmeans(m6, ~ axis_bin | l1_contrast))
em$beta <- em$emmean

pdf('pe_full_selective_betas_general_location_12_bins_persubj_rms.pdf', height = 3, width = 5)
ggplot(em, aes(axis_bin, beta, lty = l1_contrast, color = as.numeric(axis_bin), shape = l1_contrast)) + 
  geom_hline(yintercept = 0) +
  geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey80", show.legend = FALSE) + 
  #geom_linerange(aes(ymin = beta-SE, ymax = beta+SE, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey80", show.legend = FALSE) + 
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=as.numeric(axis_bin), linetype=NULL), position = position_dodge(width = .8), size = .4) + 
  geom_line(aes(x=as.numeric(axis_bin), shape=NULL, linetype=l1_contrast, color = NULL), position = position_dodge(width = .8)) + 
  geom_point(color="gray10", position = position_dodge(width = .8), size = 4.5) + 
  geom_point(position = position_dodge(width = .8), size = 3.5) + 
  
  # scale_color_viridis_d(option = "plasma", guide = 'none') +
  scale_color_gradientn(colors = pal, guide = 'none') +
  # scale_color_brewer(type = 'seq', guide = 'none') +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey65"), 
        panel.grid.minor = element_line(colour = "grey65"), 
        panel.background = element_rect(fill = 'grey60'), 
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2 , "lines"),
        legend.text = element_text(lineheight = 1.2),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=6))) +
  #scale_linetype(labels = c("Entropy\nselective", "Entropy\nfull")) + 
  #scale_shape(labels = c("Entropy\nselective", "Entropy\nfull")) + 
  #xlab("Post. <= Long axis position => Ant.") +
  xlab("Long axis position") + ylab("Response magnitude (AU)") +
  coord_cartesian(xlim=c(0.8,12.2))

dev.off()

## Revision analysis: PE trial fixed with different learning rates against SCEPTIC pe_max
## PE trial-wise versus state-wise
pe_compare_df$axis_bin <- cut(pe_compare_df$atlas_value, breaks = 12, labels = 1:12)
#pedf$isfixed <- as.numeric(pedf$l1_contrast=="pe_trial_fixed")
#pedf$issceptic <- as.numeric(pedf$l1_contrast=="pe_max")

# emm <- pairs(emmeans(m6, ~l1_contrast | axis_bin)) #, infer=TRUE, adjust="asymptotic")
# summary(emm, infer=TRUE, adjust="asymptotic")
# 
# summary(glht(m6, linfct=emm@linfct), test=adjusted("single-step"))
# 
# m7 <- lmer(beta ~ axis_bin * l1_contrast + (0 + isfixed|numid) + (0 + issceptic | numid), pedf)
# summary(m7)
# Anova(m7, '3')
# 
# anova(m6, m7)
# 
# em <- as_tibble(emmeans(m6, ~ axis_bin | l1_contrast))
# em$beta <- em$emmean

# m_lr_comparison <- lme(beta ~ axis_bin * l1_contrast,
#                        random=list(numid=pdDiag(~l1_contrast)), pe_compare_df, method="ML")

m_lr_comparison <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), pe_compare_df)
m_lr_comparison_2 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid/l1_contrast), pe_compare_df) #hack heterogeneity by model
m_lr_comparison_3 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid) + (1|l1_contrast), pe_compare_df) #hack heterogeneity by model

anova(m_lr_comparison, m_lr_comparison_2, m_lr_comparison_3)

emm <- pairs(emmeans(m_lr_comparison_2, ~l1_contrast | axis_bin)) #, infer=TRUE, adjust="asymptotic")
summary(emm, infer=TRUE)#, adjust="asymptotic")

#summary(glht(m6, linfct=emm@linfct), test=adjusted("single-step"))

which(emm@linfct[120,] != 0)

cvec_base <- rep(0, length(fixef(m_lr_comparison)))
names(cvec_base) <- names(fixef(m_lr_comparison))
model_comp <- c("l1_contrastpe_trial_fixed_p05", "l1_contrastpe_trial_fixed_p10", "l1_contrastpe_trial_fixed_p15", "l1_contrastpe_trial_fixed_p20")
cmat <- c()
#bin 1 is the reference condition
#pe_max (selective) is the reference condition
library(stringr)
for (b in 1:12) {
  for (m in model_comp) {
    this_vec <- cvec_base
    this_vec[m] <- -1
    if (b > 1) {
      this_vec[paste0("axis_bin", b, ":", m)] <- -1
    }
    cmat <- rbind(cmat, this_vec)
    rownames(cmat)[nrow(cmat)] <- paste0("bin_", b, "_pe_sel-pe_", str_sub(m,-3,-1))
  }
}

#not really a fair test since it penalizes across learning rates
#summary(glht(m_lr_comparison, linfct=cmat), test=adjusted("single-step"))

#by learning rate (comparable to above)
for (lr in c("p05", "p10", "p15", "p20")) {
  cat("LR: ", lr, "\n")
  csub <- cmat[grep(lr, rownames(cmat), value=TRUE),]
  #print(summary(glht(m_lr_comparison_2, linfct=csub), test=adjusted("single-step")))
  print(summary(glht(m_lr_comparison_2, linfct=csub), test=adjusted("Westfall")))
}

res <- pairs(emmeans(m_lr_comparison, ~  l1_contrast | axis_bin))
summary(glht(m_lr_comparison, linfct=res@linfct), test=adjusted("single-step")) #all pairs at once -- slices 1-3 survive

#for simplest comparisons to original analysis, perhaps just go pair by pair (SCEPTIC vs. target fixed LR)
m1_test <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), pe_compare_df %>% filter(l1_contrast %in% c("pe_max", "pe_trial_fixed_p20")) %>% droplevels())
res <- pairs(emmeans(m1_test, ~  l1_contrast | axis_bin))
summary(glht(m1_test, linfct=res@linfct), test=adjusted("single-step")) #all pairs at once -- slices 1-3 survive




##PE VS ENTROPY SELECTIVE
#df <- df_simult #swap in simultaneous model for a moment
df$axis_bin <- cut(df$atlas_value, breaks = 12, labels = 1:12)
df <- df %>% mutate(l1_contrast=ordered(l1_contrast, levels=c("v_entropy", "pe_max")))
df$ispe <- as.numeric(df$l1_contrast=="pe_max")
df$isentropy <- as.numeric(df$l1_contrast=="v_entropy")

m7 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), df)
summary(m7)
Anova(m7, '3')

m8 <- lmer(beta ~ axis_bin * l1_contrast + (0+ispe|numid) + (0+isentropy|numid), df)
summary(m8)
Anova(m8, '3')

summary(emmeans(m8, ~ axis_bin | l1_contrast), infer=TRUE)


anova(m7, m8)

em <- as_tibble(emmeans(m8, ~ axis_bin | l1_contrast))
em$beta <- em$emmean

pe_en_fig <- ggplot(em, aes(axis_bin, beta, color = as.numeric(axis_bin), shape = l1_contrast)) + 
  geom_hline(yintercept = 0) +
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey85", show.legend = FALSE) + 
  geom_linerange(aes(ymin = beta-SE, ymax = beta+SE, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .8, color="grey80", show.legend = FALSE) + 
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=as.numeric(axis_bin), linetype=NULL), position = position_dodge(width = .8), size = .4) + 
  geom_line(aes(x=as.numeric(axis_bin), shape=NULL, group=l1_contrast, color = NULL), position = position_dodge(width = .8), size=1) + 
  geom_point(color="gray10", position = position_dodge(width = .8), size = 4.5) + 
  geom_point(position = position_dodge(width = .8), size = 3.5) + 
  
  # scale_color_viridis_d(option = "plasma", guide = 'none') +
  scale_color_gradientn(colors = pal, guide = 'none') +
  # scale_color_brewer(type = 'seq', guide = 'none') +
  theme_gray(base_size=14) + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey65"), 
        panel.grid.minor = element_line(colour = "grey65"), 
        panel.background = element_rect(fill = 'grey60'), 
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2 , "lines"),
        legend.text = element_text(lineheight = 1.0),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=6)))  +
  scale_linetype(labels = c("Entropy\n(reversed)", "RPE")) + 
  scale_shape_manual(labels = c("Entropy\n(reversed)", "RPE"), values=c(16, 15)) + 
  #xlab("Post. <= Long axis position => Ant.") +
  xlab("Long axis position") + ylab("Response magnitude (AU)") +
  coord_cartesian(xlim=c(0.8,12.2)) + ggtitle("Entropy and RPE as\nseparate voxelwise predictors")

#ggsave('pe_entropy_selective_betas_general_location_12_bins_persubj_rms.pdf', pe_en_fig, height = 3, width = 5, useDingbats=FALSE)
#ggsave('pe_entropy_selective_betas_general_location_12_bins_zstat.pdf', pe_en_fig, height = 3, width = 5, useDingbats=FALSE)
ggsave('pe_entropy_simultaneous_selective_betas_general_location_12_bins_zstat.pdf', pe_en_fig, height = 3.5, width = 5, useDingbats=FALSE)

#cobble together multi-panel fig for supp.
#pe_en_fig_simult <- pe_en_fig
library(cowplot)
pdf("pe_entropy_simult_vs_indiv.pdf", width=11, height=4, useDingbats = FALSE)
plot_grid(pe_en_fig, pe_en_fig_simult, nrow=1, align="hv", axis="tblr") #+ theme(legend.position = "none")
dev.off()

## PE and Entropy first-half versus second half of each run
df_split$axis_bin <- cut(df_split$atlas_value, breaks = 12, labels = 1:12)
df_split <- df_split %>% mutate(
  side=if_else(atlas_name=="long_axis_l_2.3mm.nii.gz", "L", "R")
) %>% dplyr::select(-atlas_name)

df_split$ispe <- as.numeric(df_split$l1_contrast %in% c("pe_1h", "pe_2h"))
df_split$isentropy <- as.numeric(df_split$l1_contrast %in% c("v_entropy_1h", "v_entropy_2h"))

m7 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), df_split)
summary(m7)
Anova(m7, '3')

m8 <- lmer(beta ~ axis_bin * l1_contrast + (0+ispe|numid) + (0+isentropy|numid), df_split)
summary(m8)
Anova(m8, '3')

summary(emmeans(m8, ~ axis_bin | l1_contrast), infer=TRUE)

anova(m7, m8)

m7_lme <- lme(beta ~ axis_bin * l1_contrast,
              random=~1|numid, df_split, method="ML")

m7_lme_side <- lme(beta ~ axis_bin * l1_contrast + side,
              random=~1|numid, df_split, method="ML")


#allow per contrast random intercept
m9 <- lme(beta ~ axis_bin * l1_contrast,
          random=list(numid=pdDiag(~l1_contrast)), df_split, method="ML")

df_split <- df_split %>% mutate(
  effect=if_else(l1_contrast %in% c("v_entropy_1h", "v_entropy_2h"), "Entropy", "RPE"), 
  half=if_else(l1_contrast %in% c("v_entropy_1h", "pe_1h"), "Trials 1-25", "Trials 26-50")
)

# m11 <- lme(beta ~ axis_bin*effect*half, random=list(numid=pdDiag(~effect*half)), df_split)
# summary(m11)
# car::Anova(m11)

m11 <- lme(beta ~ axis_bin*effect*half, random=list(numid=pdDiag(~effect*half)), df_split)
summary(m11)
car::Anova(m11, type=3)
emmeans(m11, ~axis_bin*effect | half)

m11_1h <- lme(beta ~ axis_bin*effect, random=list(numid=pdDiag(~effect)), df_split %>% filter(half=="Trials 1-25"))
m11_2h <- lme(beta ~ axis_bin*effect, random=list(numid=pdDiag(~effect)), df_split %>% filter(half=="Trials 26-50"))
car::Anova(m11_1h, type=3)
car::Anova(m11_2h, type=3)

#what about a raw plot?
df_toplot <- df_split %>% group_by(axis_bin, effect, half) %>%
  summarize(xx=list(mean_cl_boot(beta))) %>% unnest(xx) %>% ungroup()
  #summarize(m_beta=mean(beta), se_beta=plotrix::std.error(beta)) %>% ungroup()

#Fig S2: use bootstrapped 95% CIs for more direct view of the data

ggplot(df_toplot, aes(x=axis_bin, y=y, ymin=ymin, ymax=ymax, color=effect)) + geom_pointrange() + geom_line(aes(x=as.numeric(axis_bin))) +
  facet_wrap(~half)

pe_en_fig <- ggplot(df_toplot, aes(x=axis_bin, y=y, color = as.numeric(axis_bin), shape = effect)) + 
  geom_hline(yintercept = 0) +
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey85", show.legend = FALSE) + 
  geom_linerange(aes(ymin = ymin, ymax = ymax, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .8, color="grey80", show.legend = FALSE) + 
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=as.numeric(axis_bin), linetype=NULL), position = position_dodge(width = .8), size = .4) + 
  geom_line(aes(x=as.numeric(axis_bin), shape=NULL, group=effect, color = NULL), position = position_dodge(width = .8), size=1) + 
  geom_point(color="gray10", position = position_dodge(width = .8), size = 4.5) + 
  geom_point(position = position_dodge(width = .8), size = 3.5) + 
  
  # scale_color_viridis_d(option = "plasma", guide = 'none') +
  scale_color_gradientn(colors = pal, guide = 'none') +
  # scale_color_brewer(type = 'seq', guide = 'none') +
  theme_gray(base_size=18) + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey65"), 
        panel.grid.minor = element_line(colour = "grey65"), 
        panel.background = element_rect(fill = 'grey60'), 
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2 , "lines"),
        legend.text = element_text(lineheight = 1.0),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=6)))  +
  #scale_linetype(labels = c("Entropy\n(reversed)", "RPE")) + 
  scale_linetype() + 
  scale_shape_manual(labels = c("Entropy\n(reversed)", "RPE"), values=c(16, 15)) + 
  #xlab("Post. <= Long axis position => Ant.") +
  xlab("Long axis position") + ylab("Response magnitude (AU)") +
  facet_wrap(~half) #+ ggtitle("Selective Maintenance 1st vs. 2nd half hippo betas")
#coord_cartesian(xlim=c(0.8,12.2))

#ggsave('pe_entropy_selective_betas_general_location_12_bins_persubj_rms.pdf', pe_en_fig, height = 3, width = 5, useDingbats=FALSE)
ggsave('1h_2h_pe_entropy_selective_betas_general_location_12_bins_zstat_boot95.pdf', pe_en_fig, height = 5, width = 10, useDingbats=FALSE)


#allow for separate intercepts as a function of bin, effect, and half (this runs forever --bad idea!)
# m11 <- lme(beta ~ axis_bin*effect*half, random=list(numid=pdDiag(~effect*half*axis_bin)), df_split)
# summary(m11)
# car::Anova(m11, type=3)



summary(m9)


anova(m7_lme, m7_lme_side, m9)

#add allowance for per-contrast l1 residual variance
m10 <- lme(beta ~ axis_bin * l1_contrast,
           random=list(numid=pdDiag(~l1_contrast)),
           weights=varIdent(form=~1|l1_contrast),
           df_split, method="ML")

summary(m10)
anova(m9, m10)

#em <- as_tibble(emmeans(m8, ~ axis_bin | l1_contrast))
em <- as_tibble(emmeans(m9, ~ axis_bin | l1_contrast))
em$beta <- em$emmean

em <- em %>% mutate(
  effect=if_else(l1_contrast %in% c("v_entropy_1h", "v_entropy_2h"), "Entropy", "RPE"), 
  half=if_else(l1_contrast %in% c("v_entropy_1h", "pe_1h"), "Trials 1-25", "Trials 26-50"),
  l1_contrast=ordered(l1_contrast, levels=c("v_entropy_1h", "v_entropy_2h", "pe_1h", "pe_2h")) #control over order on plot
)

car::Anova(m9)

em <- as_tibble(emmeans(m11, ~ axis_bin*effect | half))
em$beta <- em$emmean


pe_en_fig <- ggplot(em, aes(axis_bin, beta, color = as.numeric(axis_bin), shape = effect)) + 
  geom_hline(yintercept = 0) +
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey85", show.legend = FALSE) + 
  geom_linerange(aes(ymin = beta-SE, ymax = beta+SE, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .8, color="grey80", show.legend = FALSE) + 
  #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=as.numeric(axis_bin), linetype=NULL), position = position_dodge(width = .8), size = .4) + 
  geom_line(aes(x=as.numeric(axis_bin), shape=NULL, group=effect, color = NULL), position = position_dodge(width = .8), size=1) + 
  geom_point(color="gray10", position = position_dodge(width = .8), size = 4.5) + 
  geom_point(position = position_dodge(width = .8), size = 3.5) + 
  
  # scale_color_viridis_d(option = "plasma", guide = 'none') +
  scale_color_gradientn(colors = pal, guide = 'none') +
  # scale_color_brewer(type = 'seq', guide = 'none') +
  theme_gray(base_size=14) + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey65"), 
        panel.grid.minor = element_line(colour = "grey65"), 
        panel.background = element_rect(fill = 'grey60'), 
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2 , "lines"),
        legend.text = element_text(lineheight = 1.0),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=6)))  +
  #scale_linetype(labels = c("Entropy\n(reversed)", "RPE")) + 
  scale_linetype() + 
  #scale_shape_manual(labels = c("Entropy\n(reversed)", "RPE"), values=c(16, 15)) + 
  #xlab("Post. <= Long axis position => Ant.") +
  xlab("Long axis position") + ylab("Response magnitude (AU)") +
  facet_wrap(~half) #+ ggtitle("Selective Maintenance 1st vs. 2nd half hippo betas")
  #coord_cartesian(xlim=c(0.8,12.2))

#ggsave('pe_entropy_selective_betas_general_location_12_bins_persubj_rms.pdf', pe_en_fig, height = 3, width = 5, useDingbats=FALSE)
ggsave('1h_2h_pe_entropy_selective_betas_general_location_12_bins_zstat.pdf', pe_en_fig, height = 6, width = 8, useDingbats=FALSE)


##PE, PE FIXED, ENTROPY
#too busy for me
# adf$axis_bin <- cut(adf$atlas_value, breaks = 12, labels = 1:12)
# 
# m5 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), adf)
# summary(m5)
# Anova(m5, '3')
# 
# em <- as_tibble(emmeans(m5, ~ axis_bin | l1_contrast))
# em$beta <- em$emmean
# # going a little crazy on the details here
# library(wesanderson)
# pal = wes_palette("Zissou1", 12, type = "continuous")
# pdf('pe_trial_fixed_entropy_selective_betas_general_location_12_bins_persubj_rms.pdf', height = 3, width = 4.5)
# ggplot(em, aes(axis_bin, beta, lty = l1_contrast, color = as.numeric(axis_bin), shape = l1_contrast)) + 
#   geom_hline(yintercept = 0) +
#   geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey80", show.legend = FALSE) + 
#   #geom_linerange(aes(ymin = beta-SE, ymax = beta+SE, color=NULL, linetype=NULL), position = position_dodge(width = .8), size = .6, color="grey80", show.legend = FALSE) + 
#   #geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL, color=as.numeric(axis_bin), linetype=NULL), position = position_dodge(width = .8), size = .4) + 
#   geom_line(aes(x=as.numeric(axis_bin), shape=NULL, linetype=l1_contrast, color = NULL), position = position_dodge(width = .8)) + 
#   geom_point(color="black", position = position_dodge(width = .8), size = 4) + 
#   geom_point(position = position_dodge(width = .8), size = 3.0) + 
#   
#   # scale_color_viridis_d(option = "plasma", guide = 'none') +
#   scale_color_gradientn(colors = pal, guide = 'none') +
#   # scale_color_brewer(type = 'seq', guide = 'none') +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(colour = "grey55"), 
#         panel.grid.minor = element_line(colour = "grey55"), 
#         panel.background = element_rect(fill = 'grey50'), 
#         axis.text.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         legend.key.size = unit(2 , "lines"),
#         legend.text = element_text(lineheight = 1.1))  +
#   #scale_linetype(labels = c("RPE", "Entropy\n(reversed)")) + 
#   #scale_shape(labels = c("RPE", "Entropy\n(reversed)")) + 
#   #xlab("Post. <= Long axis position => Ant.") +
#   xlab("Long axis position") + ylab("") +
#   coord_cartesian(xlim=c(0.8,12.2))
# 
# dev.off()
# 
# 
# 
# 
