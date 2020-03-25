#!/bin/bash

set -ex

cobra_dir=/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/cobralab_atlases/mni_models/nifti


#liberal mask
#fslmaths mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 5.5 -uthr 12.5 -bin l_hipp_test_lib

#fslmaths mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 4.5 -uthr 5.5 -bin l_fornix
#fslmaths mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 5.5 -uthr 6.5 -bin l_fimbria

#more conservative mask would range from 6.5 -- 11.5
fslmaths $cobra_dir/mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 6.5 -uthr 11.5 -bin l_hipp_cobra_con

#fslmaths mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 6.5 -uthr 7.5 -bin l_ca1
#fslmaths mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 7.5 -uthr 8.5 -bin l_subiculum
#fslmaths mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 10.5 -uthr 11.5 -bin l_stratum #keep -- in the body and is part of the layered structure of hippo proper.. I think it is part WM, though
#fslmaths mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 11.5 -uthr 12.5 -bin l_alveus #this is a WM pathway that collects outputs to form fornix

#right conservative
fslmaths $cobra_dir/mni_icbm152_t1_tal_nlin_asym_09c_labels -thr 87.5 -uthr 92.5 -bin r_hipp_cobra_con

#shift to RPI before downsampling
#this handles the strange grid extend differences between mask and template using ANTS
3dresample -orient RPI -overwrite -prefix r_hipp_cobra_con.nii.gz -input r_hipp_cobra_con.nii.gz
3dresample -orient RPI -overwrite -prefix l_hipp_cobra_con.nii.gz -input l_hipp_cobra_con.nii.gz

#resample to 2.3 using ANTS
ResampleImageBySpacing 3 l_hipp_cobra_con.nii.gz l_hipp_cobra_con_2.3mm.nii.gz 2.3 2.3 2.3 0 0 1 #0=no smoothing, 0=0 voxel padding, 1=nn interp
ResampleImageBySpacing 3 r_hipp_cobra_con.nii.gz r_hipp_cobra_con_2.3mm.nii.gz 2.3 2.3 2.3 0 0 1 #0=no smoothing, 0=0 voxel padding, 1=nn interp

#shift back to LPI
3dresample -orient LPI -overwrite -prefix r_hipp_cobra_con_2.3mm.nii.gz -input r_hipp_cobra_con_2.3mm.nii.gz
3dresample -orient LPI -overwrite -prefix l_hipp_cobra_con_2.3mm.nii.gz -input l_hipp_cobra_con_2.3mm.nii.gz

3dresample -orient LPI -overwrite -prefix r_hipp_cobra_con.nii.gz -input r_hipp_cobra_con.nii.gz
3dresample -orient LPI -overwrite -prefix l_hipp_cobra_con.nii.gz -input l_hipp_cobra_con.nii.gz


#NB: this yields a different center:  3dinfo -header_line -same_all_grid r_hipp_cobra_con_2.3mm.nii.gz /gpfs/group/mnh5174/default/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii
#try using flirt instead

flirt -in l_hipp_cobra_con -ref $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm -applyisoxfm 2.3 -interp nearestneighbour -out l_hipp_cobra_con_2.3mm_flirt
flirt -in r_hipp_cobra_con -ref $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm -applyisoxfm 2.3 -interp nearestneighbour -out r_hipp_cobra_con_2.3mm_flirt

#make sure there are no tiny holes on the interior of the resampled mask (hurts the axis slice approach)
fslmaths l_hipp_cobra_con_2.3mm -fillh l_hipp_cobra_con_2.3mm
fslmaths r_hipp_cobra_con_2.3mm -fillh r_hipp_cobra_con_2.3mm

fslmaths l_hipp_cobra_con -fillh l_hipp_cobra_con
fslmaths r_hipp_cobra_con -fillh r_hipp_cobra_con

fslmaths l_hipp_cobra_con_2.3mm -add r_hipp_cobra_con_2.3mm -bin bilat_hipp_cobra_con_2.3mm
fslmaths l_hipp_cobra_con -add r_hipp_cobra_con -bin bilat_hipp_cobra_con

#Even though flirt does the weird RPI/LPI R-L extent transform to match the template, its NN algorithm is slightly poorer by eye
#fsleyes $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm

#4.5 -- 5.5 is l mammillary body, not part of hippo proper
#5.5 -- 6.5 is l fornix, which is in the hippocampus proper if we keep WM (output tract)


