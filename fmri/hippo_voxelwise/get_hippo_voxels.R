#read in masks
library(oro.nifti)
library(dplyr)

#harvard-oxford
#l_hippo <- readNIfTI("harvardoxford-subcortical_prob_Left_Hippocampus_2009c_thr50.nii.gz", reorient=FALSE)
#r_hippo <- readNIfTI("harvardoxford-subcortical_prob_Right_Hippocampus_2009c_thr50.nii.gz", reorient=FALSE)
#l_hippo_downsamp <- readNIfTI("harvardoxford-subcortical_prob_Left_Hippocampus_2009c_thr50_2.3mm.nii.gz", reorient=FALSE)
#r_hippo_downsamp <- readNIfTI("harvardoxford-subcortical_prob_Right_Hippocampus_2009c_thr50_2.3mm.nii.gz", reorient=FALSE)

#cobra atlas approach
l_hippo <- readNIfTI("l_hipp_cobra_con.nii.gz", reorient=FALSE)
r_hippo <- readNIfTI("r_hipp_cobra_con.nii.gz", reorient=FALSE)
l_hippo_downsamp <- readNIfTI("l_hipp_cobra_con_2.3mm.nii.gz", reorient=FALSE)
r_hippo_downsamp <- readNIfTI("r_hipp_cobra_con_2.3mm.nii.gz", reorient=FALSE)

#verified coordinate lookup against fsleyes
#note that the matrix position here is 1-based, but fsleyes is 0-based (i.e., fsleyes is one less)

lookup_coordinates <- function(map, index_value=1) {
  vox_mat <- which(map==1, arr.ind=TRUE)
  vox_mat <- cbind(vox_mat, t(apply(vox_mat, 1, function(r) { translateCoordinate(i=r, nim=map, verbose=FALSE) })))
  vox_mat <- as.data.frame(vox_mat) %>% setNames(c("i", "j", "k", "LR", "AP", "SI"))

  #look for most A + I voxel and most P + S voxel as anchors for axis
  #use summed quantiles to weight
  #1 - percent_rank(AP) is the most *anterior*

  vox_mat <- vox_mat %>% mutate(pct_p = percent_rank(AP), pct_i = percent_rank(SI),
    pct_a = 1 - pct_p, pct_s = 1 - pct_i,
    sum_pct_ai=pct_a + pct_i, sum_pct_ps=pct_p + pct_s)

  return(vox_mat)
}

#theta is radians CCW
transform_mat <- function(vox_mat, theta=0) {
  #rotate AP axis 42.51 degrees CW to reach HPC heading
  #remember: x is AP, y is SI
  vox_mat <- vox_mat %>% mutate(rot_AP = AP*cos(theta) + SI*sin(theta), rot_SI= -AP*sin(theta) + SI*cos(theta),
    long_axis=percent_rank(rot_AP))

  return(vox_mat)
}

write_long_axis <- function(map, vox_mat, out_name) {
  #write out hippocampal long axis map
  longaxis_map <- map
  longaxis_map[as.matrix(vox_mat[,c("i", "j", "k")])] <- vox_mat$long_axis
  longaxis_map@cal_min <- min(longaxis_map)
  longaxis_map@cal_min <- max(longaxis_map)
  datatype(longaxis_map) <- 16 #convert to float data type (not byte)
  longaxis_map@bitpix <- 32

  writeNIfTI(longaxis_map, out_name)
}


l_vox <- lookup_coordinates(l_hippo)
r_vox <- lookup_coordinates(r_hippo)
l_vox_downsamp <- lookup_coordinates(l_hippo_downsamp)
r_vox_downsamp <- lookup_coordinates(r_hippo_downsamp)

#spot check
r_vox %>% arrange(sum_pct_ai) %>% head(n=5)
r_vox %>% arrange(sum_pct_ps) %>% head(n=5)

#find the extremes on the 1mm R mask
top_r <- r_vox %>% arrange(sum_pct_ps) %>% head(n=1) %>% select(LR, AP, SI) %>% as.vector()
bottom_r <- r_vox %>% arrange(sum_pct_ai) %>% head(n=1) %>% select(LR, AP, SI) %>% as.vector()

#NB. I'm computing this with a sagittal view in mind. So the X axis of the equation is A-P; Y axis is S-I.
slope_r <- (top_r$SI - bottom_r$SI)/(top_r$AP - bottom_r$AP)

# To get intercept, fill in one of the points to solve for b
# SI = slope * AP
# 5 = slope * -40    ==>
intercept_r <- top_r$SI - slope_r * top_r$AP

#find the extremes on the 1mm L mask
top_l <- l_vox %>% arrange(sum_pct_ps) %>% head(n=1) %>% select(LR, AP, SI) %>% as.vector()
bottom_l <- l_vox %>% arrange(sum_pct_ai) %>% head(n=1) %>% select(LR, AP, SI) %>% as.vector()

#NB. I'm computing this with a sagittal view in mind. So the X axis of the equation is A-P; Y axis is S-I.
slope_l <- (top_l$SI - bottom_l$SI)/(top_l$AP - bottom_l$AP)

# To get intercept, fill in one of the points to solve for b
# SI = slope * AP
# 5 = slope * -40    ==>
intercept_l <- top_l$SI - slope_l * top_l$AP


#handy: http://www.webmath.com/equline1.html

#To compute 'long axisness', rotate A-P axis onto heading defined be the AI -> PS HPC vector
#Coordinate system rotation is given by:
# x' = x cos \theta + y sin \theta
# y' = -x sin \theta + y cos \theta
#
# Thus, we need to compute the rotation angle between the current A-P axis and the hippocampus vector
# Then, we rotate MNI coordinates onto this space (in 2-D, leaving X-Y alone), and use percentile_rank on the
# new "A-P" heading as long axisness'.
#
# NB. The equations above specify theta in radians, not degrees!
#
# NB. The equations above specify theta as a CCW angle. Here, we have a CW angle, need 2*pi - theta conversion

# Here is the y = mx + b for the A - P axis itself
# Flat line with slope of zero, intersecting i-s at zero
# AP = 0*SI + 0 #intersects SI at zero

#get to y = mx + b equation for HPC line

#equation for line along HPC axis
#PA = -11/12 * SI + -95/3

# And the angle between these lines is given by the difference in the arctans of their slopes
theta_r <- (atan(0) - atan(slope_r)) # Yields radians clockwise
cat("The right hippocampus is", theta_r * 180/pi, "degrees CW from the MNI AP axis\n")

theta_l <- (atan(0) - atan(slope_l)) # Yields radians clockwise
cat("The left hippocampus is", theta_l * 180/pi, "degrees CW from the MNI AP axis\n")

cat("We will take the mean of the left and right rotation calculations to get the best compromise position\n")
cat("The calculated mean rotation is: ", mean(c(theta_r, theta_l)) * 180/pi, "\n")

theta <- 2*pi - mean(c(theta_r, theta_l)) #should be CCW before coordinate transformation

#rotate AP axis 42.51 degrees CW to reach HPC heading
#remember: x is AP, y is SI
r_vox <- transform_mat(r_vox, theta)
l_vox <- transform_mat(l_vox, theta)
r_vox_downsamp <- transform_mat(r_vox_downsamp, theta)
l_vox_downsamp <- transform_mat(l_vox_downsamp, theta)

#spot checks
r_vox %>% arrange(long_axis) %>% select(LR,AP,SI,long_axis) %>% head(n=3)
r_vox %>% arrange(long_axis) %>% select(LR,AP,SI,long_axis) %>% tail(n=3)
r_vox %>% arrange(long_axis) %>% select(LR,AP,SI,long_axis) %>% filter(long_axis > .495 & long_axis < .505)
#r_vox %>% arrange(long_axis) %>% select(LR,AP,SI,long_axis) %>% filter(long_axis > .3 & long_axis < .31)
#r_vox %>% arrange(long_axis) %>% select(LR,AP,SI,long_axis) %>% filter(long_axis > .9 & long_axis < .91)

write_long_axis(r_hippo, r_vox, "long_axis_r_cobra_1mm")
write_long_axis(l_hippo, l_vox, "long_axis_l_cobra_1mm")
write_long_axis(r_hippo_downsamp, r_vox_downsamp, "long_axis_r_cobra_2.3mm")
write_long_axis(l_hippo_downsamp, l_vox_downsamp, "long_axis_l_cobra_2.3mm")

##TESTING AND VALIDATION
##whole image test to verify angle of rotation
## r_vox <- which(r_hippo_downsamp != 20, arr.ind=TRUE) #silly way to get all voxels

## r_vox <- cbind(r_vox, t(apply(r_vox, 1, function(r) { translateCoordinate(i=r, nim=r_hippo_downsamp, verbose=FALSE) })))
## r_vox <- as.data.frame(r_vox) %>% setNames(c("i", "j", "k", "LR", "AP", "SI"))

## #theta <- 45 * pi/180 #45 degrees CCW in radians for testing
## r_vox <- r_vox %>% mutate(rot_AP = AP*cos(theta) + SI*sin(theta), rot_SI= -AP*sin(theta) + SI*cos(theta),
##   long_axis=percent_rank(rot_AP))

## r_test <- r_hippo_downsamp
## r_test[as.matrix(r_vox[,c("i", "j", "k")])] <- r_vox$long_axis
## r_test@cal_min <- min(r_test)
## r_test@cal_min <- max(r_test)
## datatype(r_test) <- 16
## r_test@bitpix <- 32
## writeNIfTI(r_test, "long_axis_r_test_downsamp")
