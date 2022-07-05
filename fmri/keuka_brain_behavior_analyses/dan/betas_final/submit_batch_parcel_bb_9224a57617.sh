#!/bin/sh
module use /proj/mnhallqlab/sw/modules
module load r/4.1.2_depend
R CMD BATCH --no-save --no-restore /Users/localadmin/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/batch_run_parcel_bb_9223ec25cc.R
