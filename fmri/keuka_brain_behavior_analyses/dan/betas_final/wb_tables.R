
####

# whole-brain tables
schaefer_overlap <- read_csv(file.path(echange_dir, "schaefer_echange_parcel_overlap.csv")) %>%
  select(-retained, -prop_overlap, -nvox_parcel, -min, -nvox_retained) %>%
  mutate(max = round(max, 2), mean = round(mean, 2)) %>%
  rename(Glasser="MNI_Glasser_HCP_v1.0", Brainnetome="Brainnetome_1.0", Eickhoff_Zilles="CA_ML_18_MNI", Voxels = nvox_overlap)

library(officer)
doc <- read_docx()
doc <- body_add_table(doc, schaefer_overlap, style = "table_template")

print(doc, target = file.path(echange_dir, "echange_schaefer_table.docx"))


library(table1)
library(flextable)
library(magrittr)

# Create Table1 Object
tbl1<- table1(~Sepal.Length+Sepal.Width|Species,data=iris)
# Convert to flextable
t1flex(tbl1) %>% 
  save_as_docx(path="Iris_table1.docx")

