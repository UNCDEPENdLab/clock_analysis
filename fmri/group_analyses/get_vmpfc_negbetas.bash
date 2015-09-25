#!/bin/bash

set -e
#sub-brick 86 is m_rpe_pos_overall
files=($( find . -iname "*emoint_stats_REML+tlrc.HEAD" -type f ))

echo $files
sub_bricks=($( 3dinfo -verb ${files[@]} | grep -i 'm_rpe_pos_overall#0_Coef' | perl -pe 's/^.*At sub-brick #(\d+) .*/\1/' ))

echo ${sub_bricks[@]}


file_bricks=
for ((i=0; i < ${#files[@]}; i++)); do
    file_bricks="${file_bricks} ${files[i]}[${sub_bricks[i]}]"
done

echo $file_bricks

#vmpfc decrease is first cluster
3dROIstats -1DRformat -mask rpe_pos_ageeff_clustmask+tlrc'<1>' $file_bricks > rpe_pos_vmpfc_betas.1D
