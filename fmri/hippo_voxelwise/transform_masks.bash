#!/bin/bash
set -ex

applywarp -i harvardoxford-subcortical_prob_Left_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c --interp=spline \
	  -o harvardoxford-subcortical_prob_Left_Hippocampus_2009c

applywarp -i harvardoxford-subcortical_prob_Left_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm --interp=spline \
	  -o harvardoxford-subcortical_prob_Left_Hippocampus_2009c_2.3mm

applywarp -i harvardoxford-subcortical_prob_Right_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c --interp=spline \
	  -o harvardoxford-subcortical_prob_Right_Hippocampus_2009c

applywarp -i harvardoxford-subcortical_prob_Right_Hippocampus -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef \
	  -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm --interp=spline \
	  -o harvardoxford-subcortical_prob_Right_Hippocampus_2009c_2.3mm


fslmaths harvardoxford-subcortical_prob_Left_Hippocampus_2009c -thr 50 -bin harvardoxford-subcortical_prob_Left_Hippocampus_2009c_thr50 -odt char
fslmaths harvardoxford-subcortical_prob_Right_Hippocampus_2009c -thr 50 -bin harvardoxford-subcortical_prob_Right_Hippocampus_2009c_thr50 -odt char

fslmaths harvardoxford-subcortical_prob_Left_Hippocampus_2009c_2.3mm -thr 50 -bin harvardoxford-subcortical_prob_Left_Hippocampus_2009c_thr50_2.3mm -odt char
fslmaths harvardoxford-subcortical_prob_Right_Hippocampus_2009c_2.3mm -thr 50 -bin harvardoxford-subcortical_prob_Right_Hippocampus_2009c_thr50_2.3mm -odt char
