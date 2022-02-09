import nibabel as nib
import numpy as np
import pandas as pd

template = nib.load('Schaefer2018_200Parcels_7Networks_order.dscalar.nii')
tmp_data = template.get_fdata() #do operations here
new_data = np.zeros(tmp_data.shape)
print tmp_data.shape
#df = pyreadr.read_r('rt_encode_medusa_fmri_rt_base_parcelwise_fixed.rds')
#list = df[None]
df = pd.read_csv("fmri_brainbehavior_parcel_betas_ptfce.csv")
#roi_df = df.loc[df['term']=='fmri_beta:rt_vmax_lag']
roi_df = df.loc[df['term'] == 'rt_lag:fmri_beta']
#roi_df = roi_df.loc[roi_df['evt_time']==1]

for idx in roi_df['mask_value'].unique():
    new_data[tmp_data==int(idx)] = roi_df.loc[roi_df['mask_value'] == int(idx), 'statistic'].values[0]

new_cii = nib.cifti2.Cifti2Image(new_data, template.header)
#new_cii.to_filename('mri_beta_rt_vmax_lag.dscalar.nii')
new_cii.to_filename('mri_beta_rt_lag.dscalar.nii')
