library(oro.nifti)
library(dplyr)
library(tidyr)

## The FSL MNI-space Schaefer originals look the best and are what were released by the group. I've messed with using Fonov 2009c from the get-go, but the results are less optimal (with their normalization code)
## Thus, warp the 200-parcel FSL space to fonov using our in-lab warp coefficients and NN interpolation. This yields a good result by eye. It's a little precious on region boundaries (speckled at points, but defensibly so),
## but the centroids of the ROIs are definitely correct and a few straggler voxels will not hamper us"

## system(paste(c("applywarp -i original_masks/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm_ants",
##                "-o original_masks/Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants",
##                "-r ~/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c",
##                "-w ~/standard/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef --interp=nn"), collapse=" "))

## system(paste(c("applywarp -i original_masks/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm_ants",
##                "-o original_masks/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants",
##                "-r ~/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm",
##                "-w ~/standard/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef --interp=nn"), collapse=" "))

orig <- readNIfTI("original_masks/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants.nii.gz", reorient=FALSE)
labels <- read.table("original_masks/Schaefer2018_200Parcels_7Networks_order.txt") %>% select(1:2) %>%
  setNames(c("roinum", "network")) %>%
  mutate(hemi=if_else(grepl("_LH_", network), "L", "R"),
      network=sub("7Networks_[LR]H_", "", network)) %>%
  extract(col="network", into=c("network", "region"), regex="([^_]+)_(.*)")

ss136 <- labels %>% filter(! network %in% c("Vis", "SomMot")) #drop visual and somatomotor networks

mod <- orig
bad_mi <- which(!orig %in% unique(ss136$roinum))
mod[bad_mi] <- 0 #dump

writeNIfTI(mod, filename="masks/Schaefer_136_2.3mm")

#for parallel speed, divide parcellation into one mask per network
for (nn in unique(labels$network)) {
  thisnet <- labels %>% filter(network == nn)

  mod <- orig
  bad_mi <- which(!orig %in% unique(thisnet$roinum))
  mod[bad_mi] <- 0 #dump

  writeNIfTI(mod, filename=paste0("masks/Schaefer_", nn, "_2.3mm"))
 
}
