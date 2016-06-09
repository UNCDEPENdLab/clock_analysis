#!/bin/bash
set -ex

#alldirs=$( find /Volumes/Serena/MMClock/MR_Proc -iname "clock[0-9]" -ipath "*mni_5mm_wavelet*" -type d | grep -v bbr_exclude )

#for d in $alldirs; do
while read d; do
    d=$( dirname $d )
    set +x
    while [ $(jobs | wc -l) -ge 6 ]
    do
	sleep 10
    done
    set -x
    
    id=$( echo "$d" | perl -pe 's:.*/MR_Proc/([^/]+)/mni_5mm_wavelet/.*:\1:' )
    echo $id
    cd "${d}"
    refimg=$( find -L "/Volumes/Serena/MMClock/WPC-5640_MB/WPC5640_${id}" -iname  "*ref.hdr" | head -1 )
    #if [ ! -d bbr_noref ]; then
	sed -i .bak 's/-func_struc_dof 6/-func_struc_dof bbr/' .preproc_cmd && rm -f .preproc_cmd.bak
	sed -i .bak "1s:\$: -func_refimg $refimg:" .preproc_cmd && rm -f .preproc_cmd.bak
	mkdir bbr_noref || exists=1
	mv nfswudktm_clock* bbr_noref || nofiles=1
	mv func_to_struct* bbr_noref || nofiles=1
    #fi
    
    rm -f mc_target_brain* mc_target_mask* epiref_to_func* epiref_to_struct* nfswudktm_* fswudktm_* swudktm_* wudktm_* wktm_* template* .rescaling_complete .temporal_filtering_complete .warp_complete .prepare_fieldmap_complete .fmunwarp_complete .fieldmap_* .csf_* .wm_* nuisance_regressors.txt func_to_st* struct_to_func.mat fmap2epi_bbr.mat .func2struct_complete preprocessFunctional* .preprocessfunctional_complete subject_mask.nii.gz .func2struct_complete .motion_plots_complete .motion_censor_complete .smoothing_complete wktm_* wudktm_* swudktm_*
    preprocessFunctional -resume &
done < /Volumes/Serena/MMClock/MR_Proc/toreprocess
