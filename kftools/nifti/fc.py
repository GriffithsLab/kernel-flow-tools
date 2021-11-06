from nilearn import image, plotting, datasets
from nilearn.input_data import NiftiSpheresMasker
from nilearn.glm.first_level import FirstLevelModel, make_first_level_design_matrix

import numpy as np, matplotlib.pyplot as plt, pandas as pd, nibabel as nib
import glob
from nilearn.glm import threshold_stats_img
from nilearn import datasets
from nilearn import surface
from nilearn.datasets import load_mni152_template

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





