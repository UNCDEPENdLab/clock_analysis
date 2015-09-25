#!/bin/bash
set -e

if [ $# -lt 3 ]; then
    echo "Expect 3 inputs: <base directory to search> <output directory for stats files> \"<quoted wildcard expression identifying stats files to grab>\" <directory pattern to match>"
    exit 1
fi

basedir=$(cd "$1"; pwd)
outdir=$(cd "$2"; pwd)

#echo "nparams $#"
#*valueModel_emoint_stats_REML+tlrc*

if [ "$#" -eq 4 ]; then
    echo "find $basedir -iname $3 -type f -ipath $4"
    statfiles=$( find "$basedir" -iname "$3" -type f -ipath "$4" )
else
    echo "find $basedir -iname $3 -type f"    
    statfiles=$( find "$basedir" -iname "$3" -type f )
fi

#*glm_hrf_clock_preconvolve_valueModel_emoint_normalized_stats_REML+tlrc*

statdirname=glm_stat_files

if [ ! -d "$outdir/$statdirname" ]; then
    mkdir "$outdir/$statdirname"
fi

for s in $statfiles; do
    #subid=$( echo $s | perl -pe 's:^.*/MR_Raw/(\d+_\d+)/MBclock_recon.*$:\1:' )
    subid=$( echo $s | perl -pe "s:^$basedir/(\d+)_\d+/.*$:\1:" ) #leave out scan date for now for setup3dMEMA
    #subid=$( echo $s | perl -pe "s:^$basedir/(\d+[A-z]+)_\d+[A-z]+\d+/.*$:\1:" ) #leave out scan date for now for setup3dMEMA
    
    if [ ! -f "$outdir/$statdirname/${subid}_$(basename $s)" ]; then
	echo "cp $s $outdir/$statdirname/${subid}_$(basename $s)"
	cp "$s" "$outdir/$statdirname/${subid}_$(basename $s)"
    fi
done
