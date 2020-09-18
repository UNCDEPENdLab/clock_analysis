#!/bin/bash

#this script generates the integer-valued masks for the striatal regions from the Pauli PNAS meta-analysis

set -ex

src_dir=./original_masks
dest_dir=./masks
MRI_STDDIR=$HOME/standard

applywarp -i ${src_dir}/pauli_5cluster_striatum_FSLMNI152_2mm \
	  -o ${dest_dir}/pauli_5cluster_striatum_2.3mm \
	  -r ${MRI_STDDIR}/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm \
          -w ~/standard/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef --interp=nn

#this transforms original Choi 2012 mask to FSL MNI template. Just adjusts the grid, but doesn't move anything
#was run on longleaf. Start from the output
# mri_vol2vol --targ $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz --regheader \
# 	    --mov $G/lab_resources/parcellation/Choi_JNeurophysiol12_MNI152/Choi2012_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
# 	    --o ${src_dir}/Choi2012_7Networks_MNI152_1mm_TightMask.nii.gz

applywarp -i ${src_dir}/Choi2012_7Networks_MNI152_1mm_TightMask \
	  -o ${dest_dir}/Choi2012_7Networks_2.3mm_TightMask \
	  -r ${MRI_STDDIR}/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm \
          -w ~/standard/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef --interp=nn

applywarp -i ${src_dir}/Choi2012_7Networks_MNI152_1mm_LooseMask \
	  -o ${dest_dir}/Choi2012_7Networks_2.3mm_LooseMask \
	  -r ${MRI_STDDIR}/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm \
          -w ~/standard/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef --interp=nn

3dcalc -overwrite -LPI -a ${dest_dir}/Choi2012_7Networks_2.3mm_TightMask.nii.gz \
       -expr 'a*isnegative(x)' -prefix ${dest_dir}/l_striatum_tight_7Networks_2.3mm.nii.gz

3dcalc -overwrite -LPI -a ${dest_dir}/Choi2012_7Networks_2.3mm_TightMask.nii.gz \
       -expr 'a*ispositive(x)' -prefix ${dest_dir}/r_striatum_tight_7Networks_2.3mm.nii.gz

#this is accomplished using this little fixer script
R CMD BATCH --no-save --no-restore fix_striatum_rois.R

fslmaths ${dest_dir}/r_striatum_tight_7Networks_2.3mm -add 5 -thr 6 ${dest_dir}/r_striatum_tight_7Networks_2.3mm
fslmaths ${dest_dir}/l_striatum_tight_7Networks_2.3mm -add ${dest_dir}/r_striatum_tight_7Networks_2.3mm \
	 ${dest_dir}/bilateral_striatum_tight_7Networks_2.3mm
