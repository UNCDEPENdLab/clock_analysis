#!/bin/bash
set -ex


applywarp -i original_masks/harvardoxford-subcortical_prob_Left_Accumbens -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c

applywarp -i original_masks/harvardoxford-subcortical_prob_Left_Accumbens -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c_2.3mm

applywarp -i original_masks/harvardoxford-subcortical_prob_Right_Accumbens -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c

applywarp -i original_masks/harvardoxford-subcortical_prob_Right_Accumbens -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c_2.3mm

fslmaths masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c_thr50 -odt char

fslmaths masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c_thr50 -odt char

fslmaths masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c_2.3mm -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c_thr50_2.3mm -odt char

fslmaths masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c_2.3mm -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c_thr50_2.3mm -odt char

#thr 20 looks good by eye (not too thin)
fslmaths masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c -thr 20 -bin \
    masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c_thr20 -odt char

fslmaths masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c -thr 20 -bin \
    masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c_thr20 -odt char

fslmaths masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c_2.3mm -thr 20 -bin \
    masks/harvardoxford-subcortical_prob_Left_Accumbens_2009c_thr20_2.3mm -odt char

fslmaths masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c_2.3mm -thr 20 -bin \
    masks/harvardoxford-subcortical_prob_Right_Accumbens_2009c_thr20_2.3mm -odt char

exit 1

applywarp -i original_masks/harvardoxford-subcortical_prob_Left_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Left_Hippocampus_2009c

applywarp -i original_masks/harvardoxford-subcortical_prob_Left_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Left_Hippocampus_2009c_2.3mm

applywarp -i original_masks/harvardoxford-subcortical_prob_Right_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Right_Hippocampus_2009c

applywarp -i original_masks/harvardoxford-subcortical_prob_Right_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm --interp=spline \
	  -o masks/harvardoxford-subcortical_prob_Right_Hippocampus_2009c_2.3mm


fslmaths masks/harvardoxford-subcortical_prob_Left_Hippocampus_2009c -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Left_Hippocampus_2009c_thr50 -odt char
fslmaths masks/harvardoxford-subcortical_prob_Right_Hippocampus_2009c -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Right_Hippocampus_2009c_thr50 -odt char

fslmaths masks/harvardoxford-subcortical_prob_Left_Hippocampus_2009c_2.3mm -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Left_Hippocampus_2009c_thr50_2.3mm -odt char

fslmaths masks/harvardoxford-subcortical_prob_Right_Hippocampus_2009c_2.3mm -thr 50 -bin \
    masks/harvardoxford-subcortical_prob_Right_Hippocampus_2009c_thr50_2.3mm -odt char

#bilat mask for smooth in mask approach
fslmaths harvardoxford-subcortical_prob_Left_Hippocampus_2009c_thr50_2.3mm -add harvardoxford-subcortical_prob_Right_Hippocampus_2009c_thr50_2.3mm \
	 harvardoxford-subcortical_prob_Bilateral_Hippocampus_2009c_thr50_2.3mm
