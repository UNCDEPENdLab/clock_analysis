#!/bin/bash
set -xe

find /Volumes/Serena/MMClock/MR_Proc -iname "res4d.nii.gz" -ipath "*mni_5mm_wavelet/fsl_tc_nomeanunc/FEAT_LVL1_run*.feat/stats*" > copeList.txt

rm -f runSmoothnessAvg
while read f; do
    3dFWHMx -dset $(dirname $f)/res4d.nii.gz -mask $(dirname $f)/../mask.nii.gz -out $(dirname $f)/residualSmoothness_fwhm.txt 2> /dev/null >> runSmoothnessAvg
done < copeList.txt

Rscript -e "(colMeans(read.table(\"runSmoothnessAvg\")))"

