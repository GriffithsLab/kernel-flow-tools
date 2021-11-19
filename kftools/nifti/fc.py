from nilearn import image, plotting, datasets
from nilearn.input_data import NiftiSpheresMasker
from nilearn.glm.first_level import FirstLevelModel, make_first_level_design_matrix

import numpy as np, matplotlib.pyplot as plt, pandas as pd, nibabel as nib
import glob
from nilearn.glm import threshold_stats_img
from nilearn import datasets
from nilearn import surface
from nilearn.datasets import load_mni152_template

# The usual stuff
import os,sys,glob,numpy as np,pandas as pd
from matplotlib import pyplot as plt
from IPython.display import clear_output
from scipy import stats

# Neuroimaging stuff
from nilearn.surface import load_surf_data,load_surf_mesh,load_surface,vol_to_surf
from nilearn.datasets import fetch_surf_fsaverage,fetch_atlas_surf_destrieux,fetch_surf_nki_enhanced
from nilearn.plotting import plot_surf_stat_map
import nibabel as nib

from nilearn.connectome import ConnectivityMeasure
from nilearn.input_data import NiftiLabelsMasker





### seed coords from (Eggebrecht et al., 2014) ###
seed_coord_dict = {'vis':(-19.5, -102, -3), 
                   'aud':(-67.5, -27, 12), 
                   'mot':(-67.5, -12, 27), 
                   'DAN':(-58.5, -69, -6), 
                   'FPC':(-52.5, 24, 33), 
                   'DMN':(-43.5, 21, 51)}
seed_coord_ser = pd.Series(seed_coord_dict)

### making MNI mask ###
avg152_img = load_mni152_template()
avg152_dat = avg152_img.get_data()
avg152_mask_dat = (avg152_dat > 0).astype(float)
avg152_mask_img = nib.Nifti1Image(avg152_mask_dat, affine=avg152_img.affine)



"""
Volume-based FC functions 
=============================================================
"""


def fc_for_seeds(nifti_f,dothese=None,clip_vols=[],radius=10):

  nifti_img = nib.load(nifti_f)

  if len(clip_vols) > 0:
    clip_start,clip_end = clip_vols
    nifti_dat = nifti_img.get_data()
    nifti_dat_new = nifti_dat[:,:,:,0:clip_start].copy()
    nifti_dat_new = nifti_dat_new[:,:,:,-clip_end:].copy()
    
    nifti_img = nib.Nifti1Image(nifti_dat_new,affine=nifti_img.affine)



  n_scans = nifti_img.shape[3]
  t_r = 2
  #seed_coords = (0, -53, 26) # gonna add loop to have multiple later
  #z_per_seed = {}

  z_maps,z_maps_masked = {},{}

  
  if dothese: 
    seed_coord_dict_todo = {k: seed_coord_dict[k] for k in dothese}
  else:
    seed_coord_dict_todo = seed_coord_dict

  for seed_name,seed_coord in seed_coord_dict_todo.items():

    ##for j in range(len(seed_coord_ser)):
    ### estimate contrasts for each coords ###
    
    seed_masker = NiftiSpheresMasker([seed_coord], radius=radius)## #, standardize= True)
    seed_time_series = seed_masker.fit_transform(nifti_img)
    frametimes = np.linspace(0, t_r * (n_scans - 1), n_scans)
    design_matrix = make_first_level_design_matrix(frametimes, hrf_model='spm', add_regs=seed_time_series)
    contrast_arr = np.array([1] + [0] * (design_matrix.shape[1]-1))
    # contrasts = {'seed_based_glm': contrast_arr}

    first_level_model = FirstLevelModel(t_r=t_r, slice_time_ref=0)
    first_level_model = first_level_model.fit(run_imgs=nifti_img, design_matrices=design_matrix)
    
    z_map = first_level_model.compute_contrast(contrast_arr)
    z_map_masked,_ = threshold_stats_img(z_map,mask_img=avg152_mask_img)
    
    z_maps[seed_name] = z_map
    z_maps_masked[seed_name] = z_map_masked
    #seed[seed_coord_ser.index[j]] = z_map # z_map_masked
  
  #z_map_dict[i] = z
  #z_map_masked_dict[i] = z_per_seed

  return z_maps,z_maps_masked





def run_vol_parcbasedfc(nifti_obj, n_rois):
    '''
    get 4d nifti object and calculate functional connectivity based on Shaefer parcellations

    nifti_obj: 4D nifti object
    n_rois: number of parcellations to be used {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000}

    output: functional connectivity matrix (n_rois x n_rois np array)
    '''
    # importing parcellations map

    parcellation = datasets.fetch_atlas_schaefer_2018(n_rois=n_rois)

    atlas_fname = parcellation.maps
    labels = parcellation.labels

    # Defining mask and extracting info based on parcellation map

    masker = NiftiLabelsMasker(labels_img = atlas_fname, standardize=True, memory='nilean_cache', verbose=5 ) # might be able to tweak this to get better signal

    # Getting time series

    time_series = masker.fit_transform(nifti_obj)

    # Getting correlation matrices

    corr_measure = ConnectivityMeasure(kind='correlation')
    corr_matrix = corr_measure.fit_transform([time_series])[0]
    np.fill_diagonal(corr_matrix, 0)

    return corr_matrix




"""
Surface-based FC functions 
=============================================================
"""

def run_surf_seedbasedfc_ana(seed_name,seed_hemi,vol_dat_img=None,
                             surf_dat_lh=None,surf_dat_rh=None,do_plot=False,
                             verbose=False,radius=3,title=''):

  destrieux_atlas = fetch_atlas_surf_destrieux()
  parcellation_lh = destrieux_atlas['map_left']
  labels_lh = destrieux_atlas['labels']
  parcellation_rh = destrieux_atlas['map_right']
  labels_rh = destrieux_atlas['labels']

  fs5 = fetch_surf_fsaverage(mesh='fsaverage5')

  vtx_lh,tri_lh = load_surf_data(fs5.pial_left)
  vtx_rh,tri_rh = load_surf_data(fs5.pial_right)
  sulc_lh = load_surf_data(fs5.sulc_left)
  sulc_rh = load_surf_data(fs5.sulc_right)


  nvtcs_lh,nvtcs_rh = vtx_lh.shape[0],vtx_rh.shape[0]

  if vol_dat_img:
    surf_dat_lh = vol_to_surf(vol_dat_img,[vtx_lh,tri_lh],radius=radius)
    surf_dat_rh = vol_to_surf(vol_dat_img,[vtx_rh,tri_rh],radius=radius)

  if seed_hemi=='left':
    seed_labels = np.where(parcellation_lh == labels_lh.index(seed_name))[0]
    seed_timeseries = np.mean(surf_dat_lh[seed_labels], axis=0)
    nvtcs = nvtcs_lh
    
  elif seed_hemi=='right':
    seed_labels = np.where(parcellation_rh == labels_rh.index(seed_name))[0]
    seed_timeseries = np.mean(surf_dat_rh[seed_labels], axis=0)
    nvtcs_ = nvtcs_rh

  seed_bin_map = np.zeros(nvtcs, dtype=int)
  seed_bin_map[seed_labels] = 1

  stat_map_lh = np.zeros(nvtcs_lh)
  for i in range(nvtcs_lh): 
    stat_map_lh[i] = stats.pearsonr(seed_timeseries, surf_dat_lh[i])[0]
    
  stat_map_rh = np.zeros(nvtcs_rh)
  for i in range(nvtcs_rh): 
    stat_map_rh[i] = stats.pearsonr(seed_timeseries, surf_dat_rh[i])[0]

  stat_map_lh[np.where(np.mean(surf_dat_lh, axis=1) == 0)] = 0
  stat_map_rh[np.where(np.mean(surf_dat_rh, axis=1) == 0)] = 0

  if seed_hemi=='left':
    stat_map_lh[seed_labels] = 1
  elif seed_hemi=='right':
    stat_map_rh[seed_labels] = 1

  if verbose==False: clear_output()

    
  return stat_map_lh,stat_map_rh,surf_dat_lh,surf_dat_rh,seed_timeseries,vtx_lh,tri_lh,sulc_lh,vtx_rh,tri_rh,sulc_rh,fs5



def run_surf_parcbasedfc_ana(vol_dat_img=None,
                             surf_dat_lh=None,surf_dat_rh=None,do_plot=False,
                             verbose=False,radius=3,title=''):

  destrieux_atlas = fetch_atlas_surf_destrieux()
  parcellation_lh = destrieux_atlas['map_left']
  labels_lh = destrieux_atlas['labels']
  parcellation_rh = destrieux_atlas['map_right']
  labels_rh = destrieux_atlas['labels']

  fs5 = fetch_surf_fsaverage(mesh='fsaverage5')

  vtx_lh,tri_lh = load_surf_data(fs5.pial_left)
  vtx_rh,tri_rh = load_surf_data(fs5.pial_right)
  sulc_lh = load_surf_data(fs5.sulc_left)
  sulc_rh = load_surf_data(fs5.sulc_right)


  nvtcs_lh,nvtcs_rh = vtx_lh.shape[0],vtx_rh.shape[0]

  if vol_dat_img:
    surf_dat_lh = vol_to_surf(vol_dat_img,[vtx_lh,tri_lh],radius=radius)
    surf_dat_rh = vol_to_surf(vol_dat_img,[vtx_rh,tri_rh],radius=radius)

  unroi_lh = np.unique(destrieux_atlas.map_left)
  unroi_rh = np.unique(destrieux_atlas.map_right)

  nroi_lh = unroi_lh.shape[0]
  nroi_rh = unroi_rh.shape[0]

  unroi_lhrh = np.concatenate([unroi_lh, unroi_rh],axis=0)
  nt =  surf_dat_lh.shape[1]

  nroi_lhrh = unroi_lhrh.shape[0]

  surf_dat_lh_roi = np.zeros([nroi_lh,nt])
  surf_dat_rh_roi = np.zeros([nroi_rh,nt])

  for i_it,i in enumerate(unroi_lh): surf_dat_lh_roi[i_it,:] = surf_dat_lh[destrieux_atlas.map_left==i].mean(axis=0)
  for i_it,i in enumerate(unroi_rh): surf_dat_rh_roi[i_it,:] = surf_dat_rh[destrieux_atlas.map_left==i].mean(axis=0)

  surf_dat_lhrh_roi = np.concatenate([surf_dat_lh_roi, surf_dat_rh_roi],axis=0)

  surf_dat_lhrh_roi_fc = np.corrcoef(surf_dat_lhrh_roi)
  surf_dat_lhrh_roi_fc = np.nan_to_num(surf_dat_lhrh_roi_fc)

  coms_lh = np.array([vtx_lh[destrieux_atlas.map_left==i].mean(axis=0) for i in unroi_lh])
  coms_rh = np.array([vtx_rh[destrieux_atlas.map_right==i].mean(axis=0) for i in unroi_rh])
  coms_order_lh = np.argsort(coms_lh[:,1])
  coms_order_rh = np.argsort(coms_rh[:,1])
  coms_order_rh2 = coms_order_rh + 75
  coms_order_lhrh = np.concatenate([coms_order_lh,coms_order_rh2],axis=0)
  coms_order_lhrhv2 = np.concatenate([coms_order_lh,coms_order_lh+75],axis=0)

  surf_dat_lhrh_roi_fc_reorder = surf_dat_lhrh_roi_fc[coms_order_lhrh,:][:,coms_order_lhrh]
  surf_dat_lhrh_roi_fc_reorderv2 = surf_dat_lhrh_roi_fc[coms_order_lhrhv2,:][:,coms_order_lhrhv2]


  if verbose==False: clear_output()

    
  return (surf_dat_lhrh_roi,surf_dat_lhrh_roi_fc,surf_dat_lhrh_roi_fc_reorderv2,coms_order_lhrhv2,unroi_lhrh)



