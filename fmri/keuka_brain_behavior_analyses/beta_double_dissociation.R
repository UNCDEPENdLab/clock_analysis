# double dissociation of betas along the long axis

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
setwd("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/")
h <- read_csv('v_entropy_atlas_betas.csv.gz') 
h <- h %>% filter(l2_contrast=="overall") %>% mutate(beta = -scale(beta, center = F))
setwd("../../sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/")
pe <- read_csv('pe_max_atlas_betas.csv.gz')
pe <- pe %>% filter(l2_contrast=="overall") %>% mutate(beta = scale(beta, center = F))

df <- rbind(h, pe)

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
m2 <- lmer(beta ~ axis_bin * l1_contrast + (1|numid), df)
summary(m2)
Anova(m2, '3')
em <- as_tibble(emmeans(m2, ~ axis_bin | l1_contrast))
em$beta <- em$emmean
# going a little crazy on the details here
library(wesanderson)
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
m3 <- lmer(beta ~ location * l1_contrast + (1|numid), df)
summary(m3)
Anova(m3)
em <- as_tibble(emmeans(m3, ~ location | l1_contrast))
em$beta <- em$emmean
ggplot(em, aes(location, beta, color = l1_contrast)) + geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))


