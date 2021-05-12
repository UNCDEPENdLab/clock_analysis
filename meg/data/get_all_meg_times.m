addpath('/Users/hallquist/Downloads/');
tc_dir='/Users/hallquist/Data_Analysis/clock_analysis/meg/data/tc_orig_files/';
cd(tc_dir);
ff = dir('*.mat');
for ii=1:length(ff)
   get_meg_itis(ff(ii).name); 
end