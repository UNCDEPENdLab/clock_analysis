#!/bin/bash
set -ex
featl1dirs=$( find /storage/group/mnh5174_collab/MMClock/MR_Proc -iname "*.feat" -type d -ipath "*sceptic*" )
afniout=/storage/group/mnh5174_collab/MMClock/featl1_afnicombined

[ ! -d $afniout ] && mkdir $afniout

for f in ${featl1dirs}; do
    cd $f
    id=$( echo $f | perl -pe 's:.*/MR_Proc/([^/]+)/.*:\1:') 
    runnum=$( echo $f | perl -pe 's:.*/FEAT_LVL1_run(\d+)\.feat:\1:' )
    modelname=$( echo $f | perl -pe 's:.*/mni_5mm_wavelet/([^/]+)/FEA.*:\1:')
    feat_lvl1_to_afni.R -feat_dir $PWD

    [ ! -d ${afniout}/${modelname} ] && mkdir ${afniout}/${modelname}

    mv feat_aux+tlrc.HEAD ${id}_${modelname}_feat_aux_run${runnum}+tlrc.HEAD
    gzip feat_aux+tlrc.BRIK
    mv feat_aux+tlrc.BRIK.gz ${id}_${modelname}_feat_aux_run${runnum}+tlrc.BRIK.gz

    mv feat_stats+tlrc.HEAD ${id}_${modelname}_feat_stats_run${runnum}+tlrc.HEAD
    gzip feat_stats+tlrc.BRIK
    mv feat_stats+tlrc.BRIK.gz ${id}_${modelname}_feat_stats_run${runnum}+tlrc.BRIK.gz

    mv ${id}_${modelname}_feat_aux_run${runnum}+tlrc.HEAD \
	${id}_${modelname}_feat_aux_run${runnum}+tlrc.BRIK.gz \
	${id}_${modelname}_feat_stats_run${runnum}+tlrc.HEAD \
	${id}_${modelname}_feat_stats_run${runnum}+tlrc.BRIK.gz \
	${afniout}/${modelname}

done
