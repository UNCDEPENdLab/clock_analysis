#!/bin/bash

#this script generates the integer-valued masks for subcortical regions defined in Pauli et al., 2018
#Pauli, W. M., Nili, A. N., & Tyszka, J. M. (2018). A high-resolution probabilistic in vivo atlas of human subcortical brain nuclei. Scientific Data, 5(1).

set -ex

src_dir=./original_masks
dest_dir=./masks
MRI_STDDIR=$HOME/standard

#volumes=(1 2 3 4 5 6 7 9 10 11 16)
#labels=(Pu Ca Nac EXA GPe GPi SNc SNr PBP VTA STn)
volumes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
labels=(Pu Ca Nac EXA GPe GPi SNc RN SNr PBP VTA VeP HN HTH MN STn)

nvol=${#volumes[*]}

#there is some difficulty with voxels that overlap at a given threshold
#let's take the combined mask and max across voxels

combstr="fslmaths "
combstr_ants="fslmaths "
comb4d="fslmerge -t ${dest_dir}/pauli_combined_4d"
for (( i=0; i < $(( $nvol )); i++ )); do
    fslmaths "$src_dir/CIT168toMNI152_prob_atlas_bilat_1mm__(volume ${volumes[$i]})" \
	     -thr 0.1 -bin ${dest_dir}/pauli_${labels[$i]}_0p1

    #0=no smoothing, 0=0 voxel padding, 1=nn interp
    #ResampleImageBySpacing 3 ${dest_dir}/pauli_${labels[$i]}_0p1.nii.gz ${dest_dir}/pauli_${labels[$i]}_0p1_2.3mm.nii.gz 2.3 2.3 2.3 0 0 1 
    
    flirt -in ${dest_dir}/pauli_${labels[$i]}_0p1 -ref $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm \
	  -applyisoxfm 2.3 -interp nearestneighbour -out ${dest_dir}/pauli_${labels[$i]}_0p1_2.3mm_flirt

    #make sure there are no tiny holes on the interior of the resampled mask
    #convert to integer value matching original labeling
    fslmaths ${dest_dir}/pauli_${labels[$i]}_0p1_2.3mm_flirt -fillh -mul ${volumes[$i]} \
	     ${dest_dir}/pauli_${labels[$i]}_0p1_2.3mm_flirt

    combstr="$combstr ${dest_dir}/pauli_${labels[$i]}_0p1_2.3mm_flirt -add"
    #combstr_ants="$combstr_ants ${dest_dir}/pauli_${labels[$i]}_0p1_2.3mm -add"
    comb4d="$comb4d '$src_dir/CIT168toMNI152_prob_atlas_bilat_1mm__(volume ${volumes[$i]})'"
done

#drop last -add, then create combined output image
combstr="${combstr::-4} ${dest_dir}/pauli_combined_2.3mm_flirt"
eval $combstr

#combstr_ants="${combstr_ants::-4} ${dest_dir}/pauli_combined_2.3mm"
#eval $combstr_ants

eval $comb4d

#liberal combined mask
fslmaths ${dest_dir}/pauli_combined_4d -Tmax -thr 0.02 -bin -fillh ${dest_dir}/pauli_combined_unionmask
fslmaths ${dest_dir}/pauli_combined_4d -Tmax -mas ${dest_dir}/pauli_combined_unionmask ${dest_dir}/pauli_combined_max

#need the add 1 plus remask approach to get voxels where the max subbrik is 0 (caudate -- 0-based indexing)
fslmaths ${dest_dir}/pauli_combined_4d -Tmaxn -add 1 -mas ${dest_dir}/pauli_combined_unionmask ${dest_dir}/pauli_combined_integermask

#generate mask of DA midbrain regions
#this doesn't work -- don't feel like debugging it
#3dcalc -a "${dest_dir}/pauli_combined_maxindex.nii.gz" -expr 'a*amongst(7,10,11)' -prefix test.nii.gz

Rscript -e "library(oro.nifti);
  x <- readNIfTI('${dest_dir}/pauli_combined_integermask.nii.gz', reorient=FALSE);
  x@.Data <- round(x@.Data);
  mi <- which(x %in% c(7,10,11));
  x@.Data <- array(0L, dim=dim(x));
  x[mi] <- 1;
  writeNIfTI(x, filename='${dest_dir}/pauli_da_midbrain')
"

flirt -in ${dest_dir}/pauli_combined_integermask -ref $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm \
      -applyisoxfm 2.3 -interp nearestneighbour -out ${dest_dir}/pauli_combined_integermask_2.3mm

flirt -in ${dest_dir}/pauli_da_midbrain -ref $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm \
      -applyisoxfm 2.3 -interp nearestneighbour -out ${dest_dir}/pauli_da_midbrain_2.3mm
