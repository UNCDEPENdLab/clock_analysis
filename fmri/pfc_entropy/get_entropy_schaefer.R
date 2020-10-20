setwd("~/Data_Analysis/clock_analysis/fmri/pfc_entropy/")
library(oro.nifti)
library(tidyverse)
library(lme4)
atlas_betas <- read.csv("v_entropy_atlas_betas.csv.gz")

entropy_map <- readAFNI("group_maps/v_entropy-Intercept-Age_gfeat_stats+tlrc")
sb_labels <- unlist(strsplit(entropy_map@BRICK_LABS, "~", fixed=TRUE))
esb <- which(grepl("overall_Intercept_z", sb_labels))

schaefer <- readNIfTI("Schaefer_136_2.3mm.nii.gz", reorient=FALSE)

map <- entropy_map[,,,esb]
#binmask <- as.numeric(map[map > 3.09])
emask <- array(0, dim=dim(map))
emask[map > 3.09] <- 1

overlap <- schaefer*emask
o_df <- reshape2::melt(overlap)
o_sum <- o_df %>% group_by(value) %>% tally() %>%
  filter(n > 20 & value > 0) %>% rename(atlas_value=value, vox_clust=n)

#entropy_betas <- o_sum %>% left_join(atlas_betas)
atlas_betas <- atlas_betas %>% mutate(in_entropy_map=as.numeric(atlas_value %in% unique(o_sum$atlas_value)))

vfilt <- array(0, dim=dim(map))
vfilt[schaefer %in% unique(o_sum$atlas_value)] <- 1

entropy_filtered <- schaefer * vfilt
entropy_voxelwise_mean_df <- reshape2::melt(entropy_filtered)

# all the BA7 regions are included even if they had in-between signal (e.g. Schaefer 33)

labels <- read.table("original_masks/Schaefer2018_200Parcels_7Networks_order_manual.txt", sep="\t") %>% select(1:4) %>%
  setNames(c("roinum", "network", "allan", "include")) %>%
  mutate(hemi=if_else(grepl("_LH_", network), "L", "R"),
         network=sub("7Networks_[LR]H_", "", network)) %>%
  extract(col="network", into=c("network", "region"), regex="([^_]+)_(.*)") %>%
  #mutate(region_parcel = paste(region, roinum, sep = ':')) %>%
  mutate(region_parcel = paste(allan, roinum, sep = '_')) %>%
  filter(! network %in% c("Vis", "SomMot")) %>% #drop visual and somatomotor networks
  rename(atlas_value=roinum)

atlas_betas_lab <- as_tibble(inner_join(atlas_betas, labels, by = "atlas_value"))
atlas_betas_lab$parcel <- as.factor(atlas_betas_lab$atlas_value)

bdf <- atlas_betas_lab %>% filter(l2_contrast == "overall") # just the overall contrast
bdf <- bdf %>% mutate(beta_w = psych::winsor(beta, trim = 0.01))
ddf <- bdf %>% filter(network == 'DorsAttn')                # and just the DAN

ggplot(bdf, aes(network, beta, color = as.factor(atlas_value))) + geom_boxplot()
ggplot(bdf, aes(network, beta_w, color = as.factor(atlas_value))) + geom_boxplot()
ggplot(bdf %>% filter(network == "Default"), aes(parcel, beta_w, color = as.factor(parcel))) + geom_boxplot() + geom_hline(yintercept = 0)

ggplot(ddf, aes(region_parcel, beta_w, color = region_parcel, lty = hemi)) + geom_boxplot() + coord_flip() + geom_hline(yintercept = 0) + facet_wrap(~include)
ggplot(ddf, aes(beta_w, fill = hemi)) + geom_histogram(bins = 15) + facet_wrap(~region_parcel) + geom_vline(xintercept = 0)

summary(m1 <- lmer(beta_w ~ region_parcel + (1|id) , ddf))
em <- as_tibble(summary(emmeans::emmeans(m1, ~region_parcel), infer = T))
em$beta <- em$emmean
em <- em %>% inner_join(labels %>% select(region_parcel, hemi))
ggplot(em, aes(region_parcel, beta, ymin=lower.CL, ymax=upper.CL, size = t.ratio, color=hemi)) + geom_point() + geom_errorbar(size=1) +
 coord_flip() + facet_wrap(~t.ratio >2)

library(ggcorrplot)

ggcorrplot(ddf %>% select(numid, region_parcel, beta) %>% spread(key=region_parcel, value=beta) %>% 
             select(-numid) %>% cor(use="pairwise"),  hc.order = TRUE) +
  scale_fill_viridis_c()

ggcorrplot(ddf %>% filter(include==1) %>% select(numid, region_parcel, beta) %>% spread(key=region_parcel, value=beta) %>% 
             select(-numid) %>% cor(use="pairwise"),  hc.order = TRUE) +
  scale_fill_viridis_c()


####### examine factor structure of betas

wdf <- ddf %>% select(numid, region_parcel, beta_w) %>% pivot_wider(names_from = region_parcel, values_from = beta_w)


## first include the dud regions based on Schaefer DAN labels

fa_dud5 <- psych::fa(wdf %>% select(-numid), nfactors = 5)
fa_dud3 <- psych::fa(wdf %>% select(-numid), nfactors = 3)
fa_dud4 <- psych::fa(wdf %>% select(-numid), nfactors = 4)

print(fa_dud3$loadings, cut=0.3)
print(fa_dud4$loadings, cut=0.3)

# drop the duds
gdf <- ddf %>% filter(include==1) %>% select(numid, region_parcel, beta_w) %>% pivot_wider(names_from = region_parcel, values_from = beta_w)

to_factor <- gdf %>% select(-numid) 

fa3 <- psych::fa(to_factor, nfactors = 3, fm = 'ml')

print(fa3$loadings, cut=0.3)

fa2 <- psych::fa(to_factor , nfactors =2, fm = 'ml')
print(fa2$loadings, cut=0.3)


psych::fa.parallel(to_factor, nfactors = 5, fm = 'ml')

# bifactor model
bi_fa <- psych::omega(to_factor, title="DAN bifactor", nfactors = 3) #compare with
ff <- factor.scores(to_factor, bi_fa)

bi_fa3 <- fa(gdf %>% select(-numid),3,rotate="biquartimin")
bi_fa4 <- fa(gdf %>% select(-numid),4,rotate="biquartimin")
print(bi_fa)

#save 3 specific factor solution to data.frame.
bi_fa <- psych::fa(to_factor, nfactors = 4, rotate = "bifactor")
fscores <- factor.scores(to_factor, bi_fa)
entropy_fa <- cbind(gdf$numid, fscores$scores) %>% as.data.frame() %>%
  setNames(c("numid", "general_entropy", "med_par", "fef"))

write.csv(entropy_fa, file = "entropy_beta_bifactor_fscores.csv", row.names=FALSE)

tmpdf <- gdf %>% select(-numid)
names(tmpdf) <- substr(names(tmpdf), 1, 10)
k3=kmeans(tmpdf, centers = 2)
autoplot(k3, data=tmpdf)

