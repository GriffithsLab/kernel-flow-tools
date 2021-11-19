# -*- coding: utf-8 -*-
r"""
=================================
Nifti volumes - GLM analysis
=================================

"""

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import fetch_file
from kftools.nifti import ft_glm_ana

import numpy as np
import nibabel as nib
from nilearn.image import mean_img
from nilearn.surface import load_surf_mesh,load_surf_data,vol_to_surf
from nilearn import datasets

from nilearn.plotting import plot_glass_brain,plot_stat_map,plot_surf_stat_map
from matplotlib import pyplot as plt

# %%
# Grab the data
# --------------------------------------------------

data_dir='.'

fetch_file(data_dir=data_dir, filetype='kp-nii-hbo',
            site='pitch', task='ft', subid='sub010', sesid='ses01')

fetch_file(data_dir=data_dir, filetype='kp-nii-hbr',
            site='pitch', task='ft', subid='sub010', sesid='ses01')

fetch_file(data_dir=data_dir, filetype='kp-nii-evs',
            site='pitch', task='ft', subid='sub010', sesid='ses01')

nii_hbo_f = 'pitch_sub010_ft_ses01_1017-1706_kp-nii-hbo.nii.gz'

nii_hbr_f = 'pitch_sub010_ft_ses01_1017-1706_kp-nii-hbr.nii.gz'

nii_ev_f = 'pitch_sub010_ft_ses01_1017-1706_kp-nii-evs.tsv'

# %% 
# Analysis and viz
# --------------------------------------------------    

# %%
# Run GLM analysis

res_hbo =  ft_glm_ana(nii_hbo_f, nii_ev_f, out_fstr = '')
z_thrs,glm,evs,cntrsts,img = res_hbo
z = z_thrs['rightft_minus_rest']

# %%
# Glass brain plots

disp = plot_glass_brain(z,colorbar=True,threshold=5,black_bg=True)

# %%
# Slice view stat image plots 

disp = plot_stat_map(z,colorbar=True,threshold=5,cut_coords=[-20,-20,60])

# %%
# Project to surface 

fs5 = datasets.fetch_surf_fsaverage()

lh_dat = vol_to_surf(z,surf_mesh=fs5.pial_left)
rh_dat = vol_to_surf(z,surf_mesh=fs5.pial_right)

lhc = load_surf_data(fs5.sulc_left)
rhc = load_surf_data(fs5.sulc_right)

lhi_vtx,lhi_tri = load_surf_mesh(fs5.infl_left)
rhi_vtx,rhi_tri = load_surf_mesh(fs5.infl_right)
rhi_vtx_mod = rhi_vtx.copy()
rhi_vtx_mod[:,0] += 90
lrhi_vtx = np.concatenate([lhi_vtx, rhi_vtx_mod],axis=0)
lrhi_tri = np.concatenate([lhi_tri, rhi_tri+lhi_vtx.shape[0]])
lrh_dat = np.concatenate([lh_dat,rh_dat],axis=0)
lrhc = np.concatenate([lhc,rhc],axis=0)
lrhi_vtx_rot = np.zeros_like(lrhi_vtx)
lrhi_vtx_rot[:,0] = -lrhi_vtx[:,1]
lrhi_vtx_rot[:,1] = lrhi_vtx[:,0]
lrhi_vtx_rot[:,2] = lrhi_vtx[:,2]

fig, ax = plt.subplots(ncols=3,subplot_kw={'projection': '3d'},figsize=(8,3))#, f)

a = ax[0]
disp = plot_surf_stat_map([lrhi_vtx_rot,lrhi_tri],lrh_dat, hemi='right', view='dorsal',
                          bg_map=lrhc,threshold=5,axes=a,colorbar=False)#

a = ax[1]
disp = plot_surf_stat_map([lhi_vtx,lhi_tri],lh_dat, hemi='left', view='lateral',
                          bg_map=lhc,threshold=5,axes=a,colorbar=False)#

a = ax[2]
disp = plot_surf_stat_map([rhi_vtx,rhi_tri],rh_dat, hemi='right', view='lateral',
                          bg_map=rhc,threshold=5,axes=a,colorbar=False)#


