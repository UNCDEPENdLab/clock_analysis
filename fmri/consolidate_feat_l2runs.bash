#!/bin/bash
set -ex
featl2dirs=$( find /storage/group/mnh5174_collab/MMClock/MR_Proc -iname "*.gfeat" -type d -ipath "*sceptic*" )
afniout=/storage/group/mnh5174_collab/MMClock/featl2_afnicombined

[ ! -d $afniout ] && mkdir $afniout

for f in ${featl2dirs}; do
    cd $f
    id=$( echo $f | perl -pe 's:.*/MR_Proc/([^/]+)/.*:\1:') 
    modelname=$( echo $f | perl -pe 's:.*/mni_5mm_wavelet/([^/]+)/FEA.*:\1:')

    #make directory if needed
    [ ! -d ${afniout}/${modelname} ] && mkdir ${afniout}/${modelname}

    #copy the design matrix to the output directory, too
    cp $(dirname $PWD)/designmatrix.RData ${afniout}/${modelname}/${id}_${modelname}_designmatrix.RData

    if [ -f "${afniout}/${modelname}/${id}_${modelname}_gfeat_stats+tlrc.HEAD" ]; then
	echo "File exists: ${afniout}/${modelname}/${id}_${modelname}_gfeat_stats+tlrc.HEAD. Skipping"
	continue
    fi

    feat_lvl2_to_afni.R -gfeat_dir $PWD

    ##mv feat_aux+tlrc.HEAD ${id}_${modelname}_feat_aux_run${runnum}+tlrc.HEAD
    ##gzip feat_aux+tlrc.BRIK
    ##mv feat_aux+tlrc.BRIK.gz ${id}_${modelname}_feat_aux_run${runnum}+tlrc.BRIK.gz

    mv gfeat_stats+tlrc.HEAD ${id}_${modelname}_gfeat_stats+tlrc.HEAD
    mv gfeat_stats+tlrc.BRIK.gz ${id}_${modelname}_gfeat_stats+tlrc.BRIK.gz

    mv ${id}_${modelname}_gfeat_stats+tlrc.HEAD \
	${id}_${modelname}_gfeat_stats+tlrc.BRIK.gz \
	${afniout}/${modelname}

done
