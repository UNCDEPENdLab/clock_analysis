#!/bin/sh
module use /proj/mnhallqlab/sw/modules
module load r/4.1.2_depend
R CMD BATCH --no-save --no-restore /Users/alexdombrovski/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/batch_run_parcel_bb_4df525953424.R
