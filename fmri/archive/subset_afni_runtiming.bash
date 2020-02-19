#!/bin/bash
#runs 6 7 8 for 10711 have high movement and should be excluded
#bash-3.2$ wc -l run{1,2,3,4,5}_clock.1D
#  275 run1_clock.1D
#  277 run2_clock.1D
#  272 run3_clock.1D
#  275 run4_clock.1D
#  274 run5_clock.1D
# 1373 total

#truncate to first 1373 volumes
nvol=1373

for v in clock feedback ev rpe_pos rpe_neg mean_uncertainty rel_uncertainty; do
    for e in _  _happy_ _scram_ _fear_; do
	head -n $nvol ${v}${e}concat.1D > ${v}${e}concat_run12345.1D
    done
done

head -n $nvol ../motion_pcs.txt > ../motion_pcs_run12345.txt
head -n $nvol ../censor_intersection_concat.1D > ../censor_intersection_concat_run12345.1D
