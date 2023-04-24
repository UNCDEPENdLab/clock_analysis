#!/bin/bash

# set -ex
unsmooth_files=$( find /proj/mnhallqlab/studies/MMClock/MR_Proc -ipath "*mni_nosmooth_aroma*" -iname "fawuktm_clock[0-9].nii.gz" )

#mfile=/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/bilat_hipp_cobra_con_2.3mm.nii.gz
mfile=/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/vtavship_mask_2.3mm.nii.gz

for u in $unsmooth_files; do

    while [ $( jobs -p | wc -l ) -ge 16 ]
    do
	sleep 2
    done

    #3dBlurInMask looks for isolated voxels in the mask and removes them before smoothing
    #While principled, this also changes our mask. Instead, ask 3dBlurInMask to output the original values
    #outside the mask, then use fslmaths to mask the result. This will preserve the original (unsmoothed)
    #estimates in the 3 voxels truncated by 3dBlurInMask

    ofile=${u/.nii.gz/_vtavshipp_blur.nii.gz}
    nfile="$( dirname $ofile)/n$( basename $ofile )"
    if [ -r "$nfile" ]; then
	echo "file exists: $nfile"
	continue # skip existing files
    fi
    
    (
	echo "Blurring dataset: $u"
	3dBlurInMask -overwrite -mask $mfile \
	 	     -FWHM 5 -preserve \
	 	     -input $u -prefix $ofile >/dev/null 2>/dev/null

	fslmaths $ofile -mas $mfile $ofile #remask

	tmean=$( mktemp -u )
	fslmaths $ofile -Tmean $tmean -odt float
	fslmaths $ofile -mul 100 -div $tmean $nfile -odt float
    ) &
    
done
