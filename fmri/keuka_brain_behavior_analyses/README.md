---
title: "Code for brain-to-behavior analyses"
output: html_document
---
This folder contains most of the key inferential analyses from the paper.

Dombrovski, Luna & Hallquist, Nature Communications, 2020 \

Scripts as ordered in the paper:

# beta_double_dissociation.R
- contrast prediction error and entropy signals along the hippocampal \
long axis as presented in Fig. 2b, d

# beta_cluster_import_pca_clean.R
- import fMRI hippocampal cluster data (regression coefficients, "betas") from \
conventional whole-brain analyses \
- perform PCA on cluster betas \
- merge cluster data with behavioral data \

# beta_cluster_behavior_analyses_clean.R
- multi-level models (MLM) predicting choice based on reinforcement and hippocampal signals \
(Fig. 2e-h)
- sensitivity analyses controlling for trial, contingency, uncertainty,\
subject's performance
- both sets of analyses performed on the original (fMRI) and replication (MEG) session data

# load_medusa_data.R
- import data for Mixed-Effects DeconVolved Signal Analyses (MEDUSA)
- format these data for analyses

# medusa_event_locked_lmer_clean.R
- runs 'load_medusa_data.R'
- mixed-effects analyses of deconvolved data, including:
- RT(Vmax) - aligned ramps (Fig. 3)
- responses to reinforcement (Fig. 4 a-d)
- by running 'load_medusa_data.R' this also prepares data for \
 decoding analyses (Fig. 4 e-g)

# hipp_decon_rt_prediction.R
- with 'decode = T', main analysis for Fig. 4 e-g
- with 'u = T', exploratory analysis attempting to predict the uncertainty of the next choice \
with hippocampal activity
