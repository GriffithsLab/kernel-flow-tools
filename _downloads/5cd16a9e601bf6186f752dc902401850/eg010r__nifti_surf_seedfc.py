# -*- coding: utf-8 -*-
r"""
=================================
Nifti surface seed-based FC analyses
=================================

"""

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import fetch_file
from kftools.nifti import run_surf_seedbasedfc_ana

from nilearn.datasets import fetch_surf_nki_enhanced
from nilearn.surface import load_surf_data
import nibabel as nib

from matplotlib import pyplot as plt
from nilearn.plotting import plot_surf_stat_map


# %%
# Grab the data
# --------------------------------------------------

data_dir='.'

fetch_file(data_dir=data_dir, filetype='kp-nii-hbo',
            site='snic', task='rec', subid='sub008', sesid='ses02')
nii_hbo_f = 'snic_sub008_rec_ses02_0920-1711_kp-nii-hbo.nii.gz'

nki_dataset = fetch_surf_nki_enhanced(n_subjects=1)
nki_rsfmri_lh = load_surf_data(nki_dataset['func_left'][0])
nki_rsfmri_rh = load_surf_data(nki_dataset['func_right'][0])

# %% 
# Plotting function

def do_plot(surf_lh,surf_rh,stat_map_lh,stat_map_rh,sulc_lh,sulc_rh,title='',
            thrs=[.1,.1,.1,.1],vmaxs = [1,1,1,1]):

    fig, ax = plt.subplots(ncols=4, figsize=(18,3),subplot_kw={'projection': '3d'})

    a = ax[0]
    disp = plot_surf_stat_map([vtx_lh,tri_lh], stat_map=stat_map_lh,hemi='left',
                              view='medial', colorbar=True,bg_map=sulc_lh,bg_on_data=True,
                              cmap='cold_hot', threshold=thrs[0],vmax=vmaxs[0],axes=a)
    a = ax[1]
    disp = plot_surf_stat_map([vtx_lh,tri_lh], stat_map=stat_map_lh,hemi='left',
                              view='lateral', colorbar=True,bg_map=sulc_lh,bg_on_data=True,
                              cmap='cold_hot', threshold=thrs[1],vmax=vmaxs[1],axes=a)

    a = ax[2]
    disp = plot_surf_stat_map([vtx_rh,tri_rh], stat_map=stat_map_rh,hemi='right',
                              view='lateral', colorbar=True,bg_map=sulc_rh,bg_on_data=True,
                              cmap='cold_hot', threshold=thrs[2],vmax=vmaxs[2],axes=a)

    a = ax[3]
    disp = plot_surf_stat_map([vtx_rh,tri_rh], stat_map=stat_map_rh,hemi='right',
                              view='medial', colorbar=True,bg_map=sulc_rh,bg_on_data=True,
                              cmap='cold_hot', threshold=thrs[3],vmax=vmaxs[3],axes=a)

    fig.suptitle(title)



# %% 
# Analysis and viz
# --------------------------------------------------    

# %%
# Run FC analysis

seed_name = b'G_cingul-Post-dorsal'; 
seed_hemi = 'left'

img = nib.load(nii_hbo_f)
res_hbo = run_surf_seedbasedfc_ana(seed_name,seed_hemi,vol_dat_img=img,radius=12.0)

res_fmri = run_surf_seedbasedfc_ana(seed_name,seed_hemi,
                                    surf_dat_lh=nki_rsfmri_lh,surf_dat_rh=nki_rsfmri_rh)


# %%
# Plot

(stat_map_lh,stat_map_rh,surf_dat_lh,surf_dat_rh,_,
 vtx_lh,tri_lh,sulc_lh,vtx_rh,tri_rh,sulc_rh,_) = res_hbo
do_plot([vtx_lh,tri_lh],[vtx_rh,tri_rh],stat_map_lh,stat_map_rh,sulc_lh,sulc_rh,title='',
        thrs=[.1,.1,.1,.1],vmaxs=[.4,.4,.4,4])

(stat_map_lh,stat_map_rh,surf_dat_lh,surf_dat_rh,_,
 vtx_lh,tri_lh,sulc_lh,vtx_rh,tri_rh,sulc_rh,_) = res_fmri
do_plot([vtx_lh,tri_lh],[vtx_rh,tri_rh],stat_map_lh,stat_map_rh,sulc_lh,sulc_rh,title='',
        thrs=[.2,.2,.2,.2],vmaxs=[.8,.8,.8,8])

