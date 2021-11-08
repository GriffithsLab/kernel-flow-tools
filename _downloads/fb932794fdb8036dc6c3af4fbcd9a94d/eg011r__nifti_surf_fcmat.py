# -*- coding: utf-8 -*-
r"""
=================================
Nifti surface FC matrix analyses
=================================

"""

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import fetch_file
from kftools.nifti import run_surf_parcbasedfc_ana

import nibabel as nib
from matplotlib import pyplot as plt
import seaborn as sns

# %%
# Grab the data
# --------------------------------------------------

data_dir='.'

fetch_file(data_dir=data_dir, filetype='kp-nii-hbo',
            site='snic', task='rec', subid='sub008', sesid='ses02')
nii_hbo_f = 'snic_sub008_rec_ses02_0920-1711_kp-nii-hbo.nii.gz'


# %% 
# Analysis and viz
# --------------------------------------------------    

# %%
# Run analyses

img = nib.load(nii_hbo_f)
res = run_surf_parcbasedfc_ana(vol_dat_img=img)

surf_dat_lhrh_r,surf_dat_lhrh_roi_fc,surf_dat_lhrh_roi_fc_reorderv2,coms_order_lhrhv2,unroi_lhrh = res


# %%
# Brain plots
fig, ax = plt.subplots()
sns.heatmap(surf_dat_lhrh_roi_fc_reorderv2,cmap='cold_hot',vmin=-1,vmax=1,axes=ax)#, cmap='hot');

