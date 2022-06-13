library(tidyverse)
library(officer)

setwd("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")
betas_dir <- "~/GDrive/SCEPTIC_fMRI/wholebrain_betas"
Sys.setenv(AFNIDIR="/Users/hallquist/abin")

echange_dir <- file.path(betas_dir, "L1m-echange")

####

# whole-brain tables
schaefer_overlap <- read_csv(file.path(echange_dir, "zstat6_ptfce_Schaefer_244_final_2009c_2.3mm_overlap.csv")) %>%
  select(-prop_overlap, -nvox_parcel, -min) %>%
  mutate(max = round(max, 2), mean = round(mean, 2)) %>%
  rename(roi_num=roi_val)

# %>%
  #rename(Glasser="MNI_Glasser_HCP_v1.0", Brainnetome="Brainnetome_1.0", Eickhoff_Zilles="CA_ML_18_MNI", Voxels = nvox_overlap)

schaefer_labels <- read_csv(file.path(betas_dir, "region_labels.csv"), comment = "#") %>%
  mutate(network=recode(network, 
                        Vis="Visual",
                        SomMat = "Somatomotor",
                        DorsAttn="Dorsal Attention",
                        SalVentAttn = "Ventral Attention",
                        Cont="Frontoparietal Control"
  ))
schaefer_anat <- read_csv(file.path(betas_dir, "region_lookup.csv"))

schaefer_echange <- schaefer_overlap %>%
  inner_join(schaefer_labels, by="roi_num") %>%
  inner_join(schaefer_anat, by="roi_num") %>%
  filter(retained==TRUE) %>%
  mutate(roi_label = if_else(roi_num > 200, subregion, MNI_Glasser_HCP_v1.0)) %>%
  mutate(roi_label = sub("^\\s*(Focus point:\\s*|\\* Within \\d+ mm:\\s*)[LR]_", "", roi_label, perl=T)) %>%
  mutate(roi_label = gsub("_", " ", roi_label)) %>%
  arrange(network, hemi) %>%
  select(network,hemi,x,y,z,nvox_retained,roi_label,mean) %>%
  rename("Network"=network, "n Voxels"=nvox_retained, Hemisphere="hemi", "Label"=roi_label, "Mean z"=mean)



doc <- read_docx()
doc <- body_add_table(doc, schaefer_echange, style = "table_template")

print(doc, target = file.path(echange_dir, "echange_schaefer_table.docx"))


library(table1)
library(flextable)
library(magrittr)

# Create Table1 Object
#tbl1<- table1(~Sepal.Length+Sepal.Width|Species,data=iris)
tbl1<- table1(data=iris)
# Convert to flextable
t1flex(tbl1) %>% 
  save_as_docx(path="Iris_table1.docx")

