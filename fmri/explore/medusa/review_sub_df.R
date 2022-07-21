# examine Explore subject-level data
library(tidyverse)
library(corrplot)

setwd("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/explore_medusa/data/")

sub_df <- readRDS("explore_n146.rds")

setwd("../")

cormat <- psych::corr.test(sub_df %>% select_if(is.numeric))

pdf("explore_146_corrplot.pdf", height = 12, width = 12)
corrplot(cormat$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = cormat$p, sig.level=0.05, insig = "blank")
dev.off()
