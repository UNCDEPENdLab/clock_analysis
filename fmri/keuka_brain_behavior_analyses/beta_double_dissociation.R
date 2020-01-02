# double dissociation of betas along the long axis
library(wesanderson)
library(modelr)
library(tidyverse)
library(lme4)
library(afex)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)
library(viridis)
library(emmeans)
library(RColorBrewer)

# make long dataframe of both beta flavors
#selective maintenance

#zstat
h <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_atlas_zstats.csv.gz') %>%
  filter(l2_contrast=="overall") %>% mutate(zstat=-1*zstat) %>% dplyr::rename(beta=zstat)
#range(h$zstat)
#mutate(beta=-1*beta) 
#mutate(beta = -1*scale(beta, center = T))

#raw betas
# h <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_atlas_betas.csv.gz') %>% filter(l2_contrast=="overall") %>%
#    mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
#    group_by(ID) %>% mutate(beta = -1*scale(beta, center = F)) %>% ungroup() #per subject rms transformation

#mutate(beta = -1*scale(beta, center = F)) %>% 


# pdf("entropy_hist_by_bin_selective.pdf", width=12, height=12)
# h$axis_bin <- cut(h$atlas_value, breaks = 12, labels = 1:12)
# lattice::histogram(~zstat | axis_bin, h)
# dev.off()


#pe <- read_csv('pe_max_atlas_betas.csv.gz')
# pe <- read_csv('~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_zstats.csv.gz') %>%
#   filter(l2_contrast=="overall") %>% mutate(zstat=psych::winsor(zstat, trim=0.001)) %>% dplyr::rename(beta=zstat)
# hist(pe$zstat)
pe <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_betas.csv.gz") %>%
  filter(l2_contrast=="overall") %>%
  mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(ID) %>% mutate(beta = scale(beta, center = F)) %>% ungroup() #per subject rms transformation

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
#  group_by(ID) %>% mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = F))) %>% ungroup()
#mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = F))) %>%


#mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = F))) %>% rename(ID = id)
#mutate(l1_contrast="v_entropy_fixed", beta = -1*as.vector(scale(beta, center = T)))

pe_fixed <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_trial_fixed-preconvolve_fse_groupfixed/pe_trial_fixed/pe_trial_fixed_atlas_betas.csv.gz")
pe_fixed <- pe_fixed %>% filter(l2_contrast=="overall") %>% 
  rename(ID = id) %>%
  mutate(beta=psych::winsor(beta, trim=0.005)) %>% #1% winsorize before normalization
  group_by(ID) %>% mutate(beta = as.vector(scale(beta, center = F))) %>% ungroup()
  
hdf <- rbind(h, h_full) #%>% dplyr::rename(beta=zstat)
df <- rbind(h, pe)
adf <- rbind(df, pe_fixed)
pedf <- rbind(pe, pe_fixed)

# there is a heteroscedasticity with greater variance anteriorly
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/h_pe_betas')
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
m5 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), pedf)
summary(m5)
Anova(m5, '3')
# 100 bins
df$location <- as.factor(round(df$atlas_value, digits = 2))
m6 <- lmer(beta ~ location * l1_contrast + (1|numid), df)
summary(m6)
Anova(m6)
 em <- as_tibble(emmeans(m3, ~ location | l1_contrast))
em$beta <- em$emmean
ggplot(em, aes(location, beta, color = l1_contrast)) + geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))


setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/h_pe_betas')

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
pal = wes_palette("Zissou1", 12, type = "continuous")

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
library(multcomp)
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


##PE VS ENTROPY SELECTIVE
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
  coord_cartesian(xlim=c(0.8,12.2))

#ggsave('pe_entropy_selective_betas_general_location_12_bins_persubj_rms.pdf', pe_en_fig, height = 3, width = 5, useDingbats=FALSE)
ggsave('pe_entropy_selective_betas_general_location_12_bins_zstat.pdf', pe_en_fig, height = 3, width = 5, useDingbats=FALSE)

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
