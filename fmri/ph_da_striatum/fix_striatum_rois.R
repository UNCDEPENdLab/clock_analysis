library("oro.nifti")

##Left
x <- readNIfTI("masks/l_striatum_tight_7Networks_2.3mm.nii.gz", reorient=FALSE)

x[x==2] <- 1 #2 -> 1
highvals <- x[x > 2]
highvals <- highvals - 2 # 4 ->2; 5 -> 3; etc.
x[x > 2] <- highvals

print(table(x))

## fix range in header
x@cal_min <- min(x)
x@cal_max <- max(x)

writeNIfTI(x, "masks/l_striatum_tight_7Networks_2.3mm", verbose = TRUE)

##Right
x <- readNIfTI("masks/r_striatum_tight_7Networks_2.3mm.nii.gz", reorient=FALSE)
x[x==3] <- 0 #only 9 voxels to work with

#same renumbering as above
x[x==2] <- 1 #2 -> 1
highvals <- x[x > 2]
highvals <- highvals - 2 # 4 ->2; 5 -> 3; etc.
x[x > 2] <- highvals

print(table(x))

## fix range
x@cal_min <- min(x)
x@cal_max <- max(x)
writeNIfTI(x, "masks/r_striatum_tight_7Networks_2.3mm", verbose = TRUE)
