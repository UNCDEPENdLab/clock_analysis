#!/bin/bash
#set -x
set -e

target=$1

if [ -z ${target} ]; then
    echo "target directory not passed in as first parameter."
    echo "assuming current directory."
    target="."
fi

cd "${target}"
#allHdr=$(find $target -iname "IC*.HEAD")
allHdr=$(find "$target" -iname "ic*_onesampMixed.nii.gz")

#echo $allHdr

#FWET=5.223305 #.05
#FWET=5.643226 #.01
#FWET=6.229551 #.001

#2014 update from SPM output
#FWET=5.543850 #.001
#FWET=6.719572 #1e-06

FWET=10.0
#FWET=7.417891 #1e-08

minVoxInClusters=20

export AFNI_WHEREAMI_MAX_FIND=10

#if [ ! -d "${target}/1mm_resample" ]; then
#    mkdir "${target}/1mm_resample"
#fi

for hdr in $allHdr; do
    echo "hdr: $hdr"
    imgNum="$( echo "$hdr" | perl -pe 's/^.*ic(\d+)_onesampMixed\.nii\.gz$/\1/' )"
    #noSuffix=$( echo "$hdr" | perl -pe 's/^(.*)\.HEAD$/\1/' )

    #l/r hemi mask
    3dcalc -overwrite -LPI -a "$hdr[0]" -expr 'ispositive(x)' -prefix rmask.nii.gz
    3dcalc -overwrite -LPI -a "$hdr[0]" -expr 'isnegative(x)' -prefix lmask.nii.gz

    echo "imgNum: $imgNum"
    #3dresample -overwrite -inset ${icBase} -rmode NN -master Sci160Flip+tlrc -prefix ${imgNum}_1mm_resamp+tlrc

    tstat_subbrik=1 #second sub-brik from b, t, df, p
    3dclust -overwrite -1Dformat -noabs -nosum -1dindex $tstat_subbrik -1tindex $tstat_subbrik -1thresh ${FWET} \
	-dxyz=1 1 $minVoxInClusters "$hdr" > "ic${imgNum}_clusters.1D" #-1noneg

    3dmerge -overwrite -dxyz=1 -1clust_order 1 $minVoxInClusters \
	 -1thresh ${FWET} -1tindex $tstat_subbrik -1dindex ${tstat_subbrik} \
	-prefix "ic${imgNum}_clustMask.nii.gz" "$hdr" #-1noneg

    #look for all extrema within cluster mask
    #first 11 lines are header
    echo "#Index      Intensity    LR[mm]    PA[mm]    IS[mm]     Count  Dist[mm]" > ic${imgNum}_extrema.1D
    3dExtrema -mask_file ic${imgNum}_clustMask.nii.gz -sep_dist 25 -closure -volume -remove "$hdr[${tstat_subbrik}]" | tail -n +11 >> ic${imgNum}_extrema.1D

    #generate 10mm spheres around extrema for visualization
    tail -n +2 ic${imgNum}_extrema.1D | awk '{print $3,$4,$5,NR}' > undumpextrema

    3dUndump -overwrite -master "$hdr" -xyz -srad 6 -prefix ic${imgNum}_extrema.nii.gz undumpextrema

    #look for local extrema within each cluster
    #for each cluster, look separately in each hemisphere so that nearby l/r activations don't get merged
    nclust=$( 3dBrickStat -max ic${imgNum}_clustMask.nii.gz ) #number of clusters

    maxlocal=5
    
    #header row
    echo "#Hemisphere  Cluster    Index      Intensity    LR[mm]    PA[mm]    IS[mm]     Count  Dist[mm]" > ic${imgNum}_localExtrema.1D
    for ((i = 1; i <= $nclust; i++)); do
	3dcalc -overwrite -a ic${imgNum}_clustMask.nii.gz -expr "equals(a, $i)" -prefix tmpclustmask.nii.gz
	3dcalc -overwrite -a tmpclustmask.nii.gz -b lmask.nii.gz -expr "a*b" -prefix tmpclustmask_lhemi.nii.gz
	3dcalc -overwrite -a tmpclustmask.nii.gz -b rmask.nii.gz -expr "a*b" -prefix tmpclustmask_rhemi.nii.gz

	#[ $i -eq 13 ] && exit 1

	3dcalc -overwrite -a "$hdr[${tstat_subbrik}]" -expr 'a*ispositive(a)' -prefix tmpposact.nii.gz
	3dcalc -overwrite -a "$hdr[${tstat_subbrik}]" -expr 'a*isnegative(a)' -prefix tmpnegact.nii.gz

	#look for extrema in each hemisphere separately
	#only look within left hemi if mask includes at least 20 voxels
	#3dExtrema breaks down anyway if one or two voxels are considered...
	if [ $( 3dBrickStat -count -non-zero tmpclustmask_lhemi.nii.gz ) -ge $minVoxInClusters ]; then
	    #also need to run 3dExtrema for pos and neg activation separately since we want
	    #minima for negative activation and maxima for positive

	    for f in tmpposact.nii.gz tmpnegact.nii.gz; do
		search="-maxima"
		[ $f = tmpnegact.nii.gz ] && search="-minima"
		3dExtrema -mask_file tmpclustmask_lhemi.nii.gz -sep_dist 25 $search -closure -volume -remove \
		    "$f" | tail -n +11 | head -n $maxlocal > clustExtrema_lhemi

		nextrema=$( wc -l clustExtrema_lhemi | awk '{print $1}')
		if [ $nextrema -gt 0 ]; then
		    #perl -E "say \"$i\n\" x $nextrema" > clustCol
		    seq  -f "$i" -s '\n' $nextrema > clustCol #works on mac
		    seq  -f "L" -s '\n' $nextrema > hemiCol #works on mac
		    paste hemiCol clustCol clustExtrema_lhemi >> ic${imgNum}_localExtrema.1D
		fi
	    done
	fi

	if [ $( 3dBrickStat -count -non-zero tmpclustmask_rhemi.nii.gz ) -ge $minVoxInClusters ]; then
	    #also need to run 3dExtrema for pos and neg activation separately since we want
	    #minima for negative activation and maxima for positive

	    for f in tmpposact.nii.gz tmpnegact.nii.gz; do
		search="-maxima"
		[ $f = tmpnegact.nii.gz ] && search="-minima"
		3dExtrema -mask_file tmpclustmask_rhemi.nii.gz -sep_dist 25 $search -closure -volume -remove \
		    "$f" | tail -n +11 | head -n $maxlocal > clustExtrema_rhemi

		nextrema=$( wc -l clustExtrema_rhemi | awk '{print $1}')
		if [ $nextrema -gt 0 ]; then
		    #perl -E "say \"$i\n\" x $nextrema" > clustCol
		    seq  -f "$i" -s '\n' $nextrema > clustCol #works on mac
		    seq  -f "R" -s '\n' $nextrema > hemiCol #works on mac
		    paste hemiCol clustCol clustExtrema_rhemi >> ic${imgNum}_localExtrema.1D
		fi
	    done
	fi

    done

    #generate 10mm spheres around extrema for visualization
    tail -n +2 ic${imgNum}_localExtrema.1D | awk '{print $5, $6, $7, $2}' > undumpextrema

    3dUndump -ok_1D_text -overwrite -master "$hdr" -xyz -srad 6 -prefix ic${imgNum}_localExtrema.nii.gz undumpextrema

    rm -f hemiCol clustCol clustExtrema* tmpclustmask*.nii.gz undumpextrema rmask.nii.gz lmask.nii.gz \
	tmpposact.nii.gz tmpnegact.nii.gz

    sleep 2

    #1,2,3 sub-bricks are center of mass for x,y,z
    whereami -lpi -space MNI -coord_file "ic${imgNum}_clusters.1D[1,2,3]" -tab -max_search_radius 5 > "ic${imgNum}_roiLookup.txt"
    whereami -lpi -space MNI -ok_1D_text -coord_file "ic${imgNum}_localExtrema.1D[4,5,6]" -tab -max_search_radius 5 > "ic${imgNum}_localExtremaLookup.txt"
    whereami -lpi -space MNI -coord_file "ic${imgNum}_extrema.1D[3,4,5]" -tab -max_search_radius 5 > "ic${imgNum}_extremaLookup.txt"

done
