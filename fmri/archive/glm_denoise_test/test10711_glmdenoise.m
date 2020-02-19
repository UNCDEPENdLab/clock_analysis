%%new approach: resample data and design onto 0.2s tr grid

r1 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r1.h5','run1'));
r2 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r2.h5','run2'));
r3 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r3.h5','run3'));
r4 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r4.h5','run4'));
r5 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r5.h5','run5'));
r6 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r6.h5','run6'));
r7 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r7.h5','run7'));
r8 = single(hdf5read('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/r8.h5','run8'));

data={r1 r2 r3 r4 r5 r6 r7 r8};
clear r1 r2 r3 r4 r5 r6 r7 r8;
%tr=1.0;
stimdur=1.5;

load '/Users/michael/Data_Analysis/clock_analysis/fmri/10711_designmats_0p2s.mat';

design=struct2cell(designs);
design_zerodur=struct2cell(designs_zerodur);

tr = 0.2; %upsampled 5x
%[results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,[],[],[],'example1figures_tr0p2');

%optimize HRF, but just use event onsets, not durations
[results_onsetonly, denoiseddata_onsetonly] = GLMdenoisedata(design_zerodur,data,stimdur,tr,'optimize',[],[],'fig10711_hrfonset');

%assume hrf of 0.2s duration and an 0.2s TR.
%get error re: singular
hrf = getcanonicalhrf(0.2,0.2);
[results_assume,denoiseddata_assume] = GLMdenoisedata(design,data,stimdur,tr,'assume',hrf,[],'example1figures_assume');

%fir model -- just use onsets
[results_fir,denoiseddata_fir] = GLMdenoisedata(design_zerodur,data,stimdur,tr,'fir',60,struct('numboots',0),'fig10711_fir');

[nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock1/fswudktm_clock1_5_3mm.nii');

nii_stat.hdr=nii.hdr;
nii_stat.hdr.dime.dim = [4 64 77 64 size(results_fir.modelmd,5) 1 1 1];
nii_stat.img = squeeze(results_fir.modelmd(:,:,:,1,:));
save_nii(nii_stat, 'clock_onset_fir')
nii_stat.img = squeeze(results_fir.modelmd(:,:,:,2,:));
save_nii(nii_stat, 'feedback_onset_fir')

save('results_fir.mat', 'results_fir');



%%NEW APPROACH: USE EXISTING CONVOLVED DESIGN MATRICES
%r1 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock1/fswudktm_clock1_5_drop6_trunc281.nii.gz');
%r2 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock2/fswudktm_clock2_5_drop6_trunc283.nii.gz');
%r3 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock3/fswudktm_clock3_5_drop6_trunc278.nii.gz');
%r4 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock4/fswudktm_clock4_5_drop6_trunc281.nii.gz');
%r5 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock5/fswudktm_clock5_5_drop6_trunc280.nii.gz');
%r6 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock6/fswudktm_clock6_5_drop6_trunc280.nii.gz');
%r7 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock7/fswudktm_clock7_5_drop6_trunc288.nii.gz');
%r8 = load_untouch_nii('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock8/fswudktm_clock8_5_drop6_trunc275.nii.gz');

%save('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/fswudktm_runs.mat', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8');

%load '/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test/fswudktm_runs.mat'
load '/Users/michael/10711_glmdenoise/fswudktm_runs.mat'
cd '/Users/michael/10711_glmdenoise';

data={r1.img r2.img r3.img r4.img r5.img r6.img r7.img r8.img};
clear r1 r2 r3 r4 r5 r6 r7 r8;
tr=1.0;
stimdur=1.0;

%load in convolved design matrix
load '/Users/michael/Data_Analysis/clock_analysis/fmri/10711_designmats_convolved.mat';
design=struct2cell(designs)'; %make designs 1 x 8 to match data

%take out mean uncertainty for first pass
for i=1:length(design)
    d = design{i};
    d(:,find(ismember(regnames,'mean_uncertainty'))) =[];
    %zscore regressors to avoid terrible scaling problems in matrix inversion for computing param estimates
    d=zscore(d);
    design{i} = d;
end

regnames=regnames(~strcmp(regnames, 'mean_uncertainty'));
fprintf('%.4f ', std(design{2}));

[results_conv,denoiseddata_conv] = GLMdenoisedata(design,data,stimdur,tr,'assume',1,[],'fig10711_preconv');

save('results_conv.mat', 'results_conv');



%%
%write out glmdenoise results


[nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock1/fswudktm_1vol.nii');

%merge beta stats, se, and t-stat into single image
allstat=zeros([size(results_conv.modelmd{2}(:,:,:,1)), size(results_conv.modelmd{2}, 4)*3]);

el=1;
for i=1:size(results_conv.modelmd{2}, 4)
   allstat(:,:,:,el) =  results_conv.modelmd{2}(:,:,:,i);
   allstat(:,:,:,el+1) =  results_conv.modelse{2}(:,:,:,i);
   allstat(:,:,:,el+2) =  results_conv.modelmd{2}(:,:,:,i)./results_conv.modelse{2}(:,:,:,i);
   el=el+3;
end

nii_stat.hdr=nii.hdr;
nii_stat.hdr.dime.dim = [4 84 100 84 size(results_conv.modelmd{2}, 4)*3 1 1 1];
nii_stat.img = allstat;
save_nii(nii_stat, 'clock_conv_allbeta')


%%old style
r1clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run1_clock_FSL3col.txt','\t');
r2clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run2_clock_FSL3col.txt','\t');
r3clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run3_clock_FSL3col.txt','\t');
r4clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run4_clock_FSL3col.txt','\t');
r5clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run5_clock_FSL3col.txt','\t');
r6clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run6_clock_FSL3col.txt','\t');
r7clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run7_clock_FSL3col.txt','\t');
r8clock=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run8_clock_FSL3col.txt','\t');

r1feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run1_feedback_FSL3col.txt','\t');
r2feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run2_feedback_FSL3col.txt','\t');
r3feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run3_feedback_FSL3col.txt','\t');
r4feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run4_feedback_FSL3col.txt','\t');
r5feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run5_feedback_FSL3col.txt','\t');
r6feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run6_feedback_FSL3col.txt','\t');
r7feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run7_feedback_FSL3col.txt','\t');
r8feedback=dlmread('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run8_feedback_FSL3col.txt','\t');

design={{r1clock(:,1) r1feedback(:,1)} ...
{r2clock(:,1) r2feedback(:,1)} ...
{r3clock(:,1) r3feedback(:,1)} ...
{r4clock(:,1) r4feedback(:,1)} ...
{r5clock(:,1) r5feedback(:,1)} ...
{r6clock(:,1) r6feedback(:,1)} ...
{r7clock(:,1) r7feedback(:,1)} ...
{r8clock(:,1) r8feedback(:,1)}};

[results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,[],[],[],'example1figures');