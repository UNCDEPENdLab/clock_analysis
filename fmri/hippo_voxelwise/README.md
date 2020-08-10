## hippo_voxelwise

This folder contains the scripts and outputs for key outputs and analyses of hippocampal signals, particularly those pertaining to the
analysis of time courses in the hippocampus (e.g., Figures 3 and 4 in Dombrovski, Luna, & Hallquist).

## core scripts

### `decon_alignment.R`

This script takes deconvolved time series produced by `extract_sceptic_atlas_deconvolved.R` and event-locks them to particular events in the experiment. This is achieved by the function `event_lock_decon`, which uses linear interpolation (but not resampling) to lock BOLD activity estimates to particular events in the design such as feedback onset. The script produces a set of files with the suffix `*decon_locked.csv.gz`, which contain time course estimates for each subject and trial (experiment trials nested within subject) for a given event. For example, the 'feedback_long' files contain time-locked activity estimates for each trial beginning one second before feedback onset and ending 10 seconds after feedback onset.

### `get_hippo_voxels.R`

This script produces the bilateral long axis gradient along the hippocampus using a priori mask files provided as inputs. More specifically, it reads in the files specified by `l_hippo` and `r_hippo`, computes the coordinates of all voxels in MNI space, then tries to find the extrema along the long axis. To do so, it calculates the percentile rank in the A->P and I->S coordinate planes, averages these percentiles, then finds the top 10 and bottom 10 voxels. The top voxels are the most posterior and superior the bottom voxels are the most anterior and inferior. The script then finds the centroids of these clusters and computes the slope of the line that connects the centroids in each hemispheric mask, then computes the average of these slopes for a bilateral gradient.

After finding the slope of the line along the hippocampal long axis, it computes a linear gradient (0..1) along the voxels of the mask, then computes the rotation of the coordinate plane from a level line to the hippocampal long axis, effectively interpolating the linear gradient onto the long axis. This results in a 0..1 value in each voxel of the hippocampal mask where voxels close to 0 are in the posterior tail and voxels close to 1 are in the anterior head. These masks are saved to NIfTI images and used in subsequent processing (e.g., `extract_sceptic_atlas_deconvolved.R`) to tag voxelwise statistics with their position along the long axis.

## other scripts

### `create_hippo_masks_cobra.bash`

This script creates the gray-matter-weighted mask of the hippocampus at 2.3mm based on the COBRA atlas from Winterburn and colleagues: http://cobralab.ca/atlases/Hippocampus-subfields/. These files were used in the more conservative analyses of the hippocampal gradient reported in the supplement of Dombrovski, Luna, and Hallquist.

### `smooth_hipp_in_mask.bash`

This script applies smoothing only within the hippocampal mask to otherwise unsmoothed preprocessed fMRI data. It yields whole-brain images where the hippocampus (outlined by the specified mask file) is smoothed within mask by AFNI 3dBlurInMask, while other voxels are unsmoothed. These data were used in several analyes of Dombrovsk, Luna, and Hallquist to mitigate the concern that activity from voxels outside of the hippocampus was being blurred/blended into hippocampal estimates.

### `transform_masks.bash`

This script thresholds and transforms Harvard-Oxford hippocampal masks into the 2.3mm MNI space used for all core analyses.

## subfolders

### hippo_brms

Not currently reported in a publication. These are Bayesian multi-level VAR analyses of deconvolved BOLD activity along the long axis of the hippocampus. The goal was to look at intrahippocampal activity within a VAR framework.

### hippo_dcm

Not currently reported in a publication. These are in-progress files for a Dynamic Causal Modeling approach to hippocampal-prefrontal connectivity.

### masks

These contain processed mask files used in the pipeline. Two files in particular are instrumental in most key analyses: `long_axis_l_2.3mm.nii.gz` and `long_axis_r_2.3mm.nii.gz`. These are the primary outputs of `get_hippo_voxels.R` and contain transformed Harvard-Oxford masks the hippocampus with the voxel values indicating position along the long axis on the 0..1 scale.

### original_masks

These contain unprocessed original masks for the hippocampus, vmPFC, and other cortical regions used for extracting voxelwise statistics from a priori ROIs.
