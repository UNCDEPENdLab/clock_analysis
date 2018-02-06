#!/bin/bash
set -xe

find /gpfs/group/mnh5174/default/MMClock/MR_Proc -iname "res4d.nii.gz" -ipath "*mni_5mm_aroma/sceptic_vchosen_ventropy_dauc_pemax_preconvolve/FEAT_LVL1_run*.feat/stats*" > copeList.txt

rm -f runSmoothnessAvg
while read f; do
    while [ $(jobs | wc -l) -gt 20 ]; do
	sleep 2;
    done; 

    (
	3dFWHMx -dset $(dirname $f)/res4d.nii.gz -mask $(dirname $f)/../mask.nii.gz -out $(dirname $f)/residualSmoothness_fwhm.txt -acf $(dirname $f)/residualSmoothness_acf.txt 2> /dev/null > $(dirname $f)/residualSmoothness_fwhm_stdout
	head -n 1 $(dirname $f)/residualSmoothness_fwhm_stdout >> runSmoothnessAvg_FWHM
	tail -n 1 $(dirname $f)/residualSmoothness_fwhm_stdout >> runSmoothnessAvg_ACF
    ) &
	
done < copeList.txt
wait

Rscript -e "(colMeans(read.table(\"runSmoothnessAvg_ACF\")))"

