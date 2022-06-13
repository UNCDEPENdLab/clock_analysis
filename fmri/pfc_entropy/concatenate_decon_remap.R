# concatenate decons

filedir <- "/proj/mnhallqlab/users/michael/sceptic_decon"
subid <- "11316"
runs <- 1:8
# files <- list.files(path = filedir, pattern = paste0("sub", subid, "_run[1-8].*_deconvolved.csv.gz"), recursive = TRUE)
files <- list.files(path = filedir, pattern = paste0("sub", subid, "_run1_.*_deconvolved.csv.gz"), recursive = TRUE, full.names=TRUE)
files <- files[grepl("Schaefer.*", files)]

library(data.table)
all_f <- rbindlist(lapply(files, fread))

# check that no xyz is duplicated
t1 <- all_f[time==2,]
t1[,xx:=paste0(x,y,z)]
aa <- duplicated(t1$xx)
table(aa)

schaefer_400 <- readNIfTI("/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/original_masks/Schaefer2018_400Parcels_7Networks_order_fonov_2.3mm_ants.nii.gz", reorient = FALSE)
schaefer_200 <- readNIfTI("/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/original_masks/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants.nii.gz", reorient = FALSE)


# get indices of mask within matrix (ijk)
zero_thresh <- 1e-4 # for binarizing/indexing
a_indices_400 <- which(schaefer_400 > zero_thresh, arr.ind = TRUE)
a_indices_200 <- which(schaefer_200 > zero_thresh, arr.ind = TRUE)

# sanity check
identical(a_indices_400, a_indices_200)

# look up spatial coordinates of voxels in atlas (xyz)
a_coordinates <- t(apply(a_indices_400, 1, function(r) {
  translateCoordinate(i = r, nim = schaefer_400, verbose = FALSE)
})) %>% round(1)

# sanity check
a_coordinates_200 <- t(apply(a_indices_200, 1, function(r) {
  translateCoordinate(i = r, nim = schaefer_200, verbose = FALSE)
})) %>% round(1)

identical(a_coordinates_200, a_coordinates)

colnames(a_coordinates) <- c("x", "y", "z")

a_coordinates <- cbind(a_coordinates, roi_400=schaefer_400[a_indices_400], roi_200=schaefer_200[a_indices_200]) %>% data.frame()
setDT(a_coordinates)

# map between old and new ROIs
schaefer_dir <- "/proj/mnhallqlab/lab_resources/parcellation/schaefer_wb_parcellation"
schaefer_7_400 <- read.csv(file.path(schaefer_dir, "labels/Schaefer2018_400Parcels_7Networks_order.csv")) %>%
  mutate(network = factor(network), net_num = as.numeric(network)) %>%
  rename(network400 = network, net_num400 = net_num, roi_400=roi_num, hemi400=hemi, subregion400=subregion)

schaefer_7_200 <- read.csv(file.path(schaefer_dir, "labels/Schaefer2018_200Parcels_7Networks_order.csv")) %>%
  mutate(network = factor(network), net_num = as.numeric(network)) %>%
  rename(network200 = network, net_num200 = net_num, roi_200 = roi_num, hemi200 = hemi, subregion200 = subregion)

a_df <- merge(a_coordinates, schaefer_7_400, by = "roi_400") %>%
  merge(schaefer_7_200, by="roi_200")

dan_nodes <- a_df %>% filter(network200 == "DorsAttn" | network400 == "DorsAttn")

xtabs(~network400 + network200, dan_nodes)
xtabs(~ network400, dan_nodes %>% filter(network200=="DorsAttn"))
xtabs(~ network200, dan_nodes %>% filter(network400=="DorsAttn"))

dan_nodes %>% count(roi_200, network200, network400)

dan_nodes %>% count(network400)
dan_nodes %>% count(network200)

# voxelwise roi mapping
dan_nodes %>% count(roi_200, roi_400)
write.csv(dan_200_400_match, file = "dan_200_400_voxel_overlap.csv", row.names = F)


#a_coordinates[, vloc := paste0(x, y, z)]
setkeyv(a_coordinates, c("x", "y", "z"))
setkeyv(all_f, c("x", "y", "z"))
refresh <- merge(all_f, a_coordinates)
refresh[, atlas_name := "Schaefer_400_remap"]
refresh[, atlas_value := NULL]

all.equal(refresh$atlas_value, refresh$roi_200)

# slow and unnecessary?
#all_f[, vloc := paste0(x, y, z)]

# scale up to all files
all_schaefer <- list.dirs("/proj/mnhallqlab/users/michael/sceptic_decon", recursive=FALSE) %>% grep("Schaefer", x=., value=TRUE)
all_files <- do.call(c, lapply(all_schaefer, function(x) {
  list.files(path = file.path(x, "deconvolved"), pattern = "sub.*\\.csv.gz", full.names = TRUE)
}))

all_df <- data.frame(
  file = all_files,
  sub = sub(".*/deconvolved/sub(\\d+).*", "\\1", all_files, perl = TRUE),
  run = sub(".*/deconvolved/sub\\d+_run(\\d+).*", "\\1", all_files, perl = TRUE)
) %>%
  group_by(sub, run) %>%
  group_split()

library(fmri.pipeline)

lapply(all_df, function(sub_data) {

  r_batch <- R_batch_job$new(
    job_name="remap400",
    n_cpus = 1, mem_per_cpu = "16g",
    wall_time = "00:10:00", scheduler = "slurm",
    input_objects = fmri.pipeline:::named_list(a_coordinates, sub_data),
    r_packages = c("fmri.pipeline", "data.table"),
    r_code = expression(
      f <- rbindlist(lapply(sub_data$file, fread)),
      setkeyv(a_coordinates, c('x', 'y', 'z')),
      setkeyv(f, c('x', 'y', 'z')),
      f <- merge(f, a_coordinates),
      f[, atlas_name := 'Schaefer_400_remap'],
      f[, atlas_value := NULL],
      out_file <- file.path('/proj/mnhallqlab/users/michael/sceptic_decon/Schaefer_400_remap/deconvolved', paste0('sub', sub_data$sub[1L], '_run', sub_data$run[1L], '_Schaefer_400_remap_2.3mm_deconvolved.csv.gz')),
      fwrite(f, file=out_file)
    )
  )

  r_batch$submit()
  Sys.sleep(0.3) # avoid fast job submission

})

