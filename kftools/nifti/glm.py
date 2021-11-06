

# Usual suspects
import os,sys,glob,numpy as np,pandas as pd
from matplotlib import pyplot as plt
from IPython.display import Image,display

# Nipy neuroimaging library stuff
import nibabel as nib

from nilearn.plotting import (plot_stat_map,plot_glass_brain,
                              plot_surf_stat_map,plot_carpet,
                              plot_design_matrix,plot_contrast_matrix)
from nilearn.image import index_img,mean_img
from nilearn import masking
from nilearn.surface import vol_to_surf
from nilearn import datasets
from nilearn.glm.first_level import FirstLevelModel
from nilearn.glm import threshold_stats_img
from nilearn.reporting import get_clusters_table


def ft_glm_ana(nii_f, ev_f, out_fstr = 'ft_glm',make_figs=False):
  
  # Load data
  img = nib.load(nii_f)
  toffset = img.header.get('toffset')

  # fsaverage5 for surface plotting
  fs5 = datasets.fetch_surf_fsaverage()

  
  # Get events
  evs = pd.read_csv(ev_f,sep='\t')
  evs.timestamp -= toffset  # adjust for toffsets
  evs = evs.set_index('timestamp').loc[0:].reset_index()
  evs['trial_type'] = 'rest'
  evs['trial_type'][evs.block_type=='left'] = 'left'
  evs['trial_type'][evs.block_type=='right'] = 'right'
  evs = evs.loc[(evs.event=='start_trial') | (evs.event=='start_rest')]
  evs = evs[['timestamp', 'duration', 'trial_type']]
  evs.columns = ['onset', 'duration', 'trial_type']

  # GLM
  glm = FirstLevelModel(t_r=1,noise_model='ar1',standardize=True,#False,
                                      hrf_model='spm',drift_model='cosine',
                                      high_pass=.01)
  glm = glm.fit(img,events=evs)


  # Contrasts
  N = glm.design_matrices_[0].shape[1]; 
  tmp = np.zeros([N])
  conditions = {k: tmp.copy() for k in ['leftft', 'rest', 'rightft']}
  conditions['leftft'][0] = 1; conditions['rest'][1] = 1; conditions['rightft'][2] = 1
  leftft_minus_rest = conditions['leftft'] - conditions['rest']
  rightft_minus_rest = conditions['rightft'] - conditions['rest']
  leftft_minus_rightft = conditions['leftft'] - conditions['rightft']
  rightft_minus_leftft = conditions['rightft'] - conditions['leftft']
  ft_minus_rest = (conditions['leftft'] +  conditions['rightft'] ) - 2 * (conditions['rest'])
  rest_minus_ft = -ft_minus_rest

  cntrsts = \
  {'leftft_minus_rest': [leftft_minus_rest,glm.compute_contrast(leftft_minus_rest,output_type='all')],
  'rightft_minus_rest': [rightft_minus_rest,glm.compute_contrast(rightft_minus_rest,output_type='all')],
  'leftft_minus_rightft': [leftft_minus_rightft,glm.compute_contrast(leftft_minus_rightft,output_type='all')],
  'rightft_minus_leftft': [rightft_minus_leftft,glm.compute_contrast(rightft_minus_leftft,output_type='all')],
  'ft_minus_rest': [ft_minus_rest,glm.compute_contrast(ft_minus_rest,output_type='all')],
  'rest_minus_ft': [rest_minus_ft,glm.compute_contrast(rest_minus_ft,output_type='all')]}

  # Thresholding
  z_thrs = {}
  for k,v in cntrsts.items():
    z = v[1]['z_score']
    z_thr, threshold = threshold_stats_img(z,alpha=.05,two_sided=False,
                                           height_control='fpr', cluster_threshold=0)
    z_thrs[k] = z_thr


  if make_figs:

    # GB plot

    fig, ax = plt.subplots(ncols=3,nrows=1, figsize=(12,3))

    a = ax[0]
    z = z_thrs['leftft_minus_rightft']
    disp = plot_glass_brain(z,axes=a,colorbar=True);
    a.set_title('leftft > rightft');

    a = ax[1]
    z = z_thrs['rightft_minus_leftft']
    disp = plot_glass_brain(z,axes=a,colorbar=True);
    a.set_title('rightft > leftft');

    a = ax[2]
    z = z_thrs['ft_minus_rest']
    disp = plot_glass_brain(z,axes=a,colorbar=True);
    a.set_title('ft > rest');

    outf = out_fstr + '_glassbrain.png'
    plt.savefig(outf)
    plt.close()


    # Slice plot

    fig, ax = plt.subplots(ncols=3,nrows=2, figsize=(23,3))

    z = z_thrs['leftft_minus_rightft']
    disp = plot_stat_map(z,black_bg=False,cut_coords = [4,6,8,10,12,14],
                        display_mode='y',axes=ax[0][0],annotate=True)
    disp = plot_stat_map(z,black_bg=False,cut_coords = [54,52,50,48,46,44],
                        display_mode='y',axes=ax[1][0],annotate=True)
    
    z = z_thrs['rightft_minus_leftft']
    disp = plot_stat_map(z,black_bg=False,cut_coords = [4,6,8,10,12,14],
                        display_mode='y',axes=ax[0][1],annotate=True)
    disp = plot_stat_map(z,black_bg=False,cut_coords = [54,52,50,48,46,44],
                        display_mode='y',axes=ax[1][1],annotate=True)

    z = z_thrs['ft_minus_rest']
    disp = plot_stat_map(z,black_bg=False,cut_coords = [4,6,8,10,12,14],
                        display_mode='y',axes=ax[0][2],annotate=True)
    disp = plot_stat_map(z,black_bg=False,cut_coords = [54,52,50,48,46,44],
                        display_mode='y',axes=ax[1][2],annotate=True)
    
    outf = out_fstr + '_slices.png'
    plt.savefig(outf)
    plt.close()

    # Surf plot


    fig, ax = plt.subplots(ncols=6,nrows=2, figsize=(30,7), subplot_kw={'projection': '3d'})


    z =  z_thrs['leftft_minus_rightft']
    dat_onsurf_lh = vol_to_surf(z, fs5.pial_left,radius=3)
    dat_onsurf_rh = vol_to_surf(z, fs5.pial_right,radius=3)

    a = ax[0][0]
    disp = plot_surf_stat_map(fs5.infl_left, dat_onsurf_lh, hemi='left', view='lateral', 
                              bg_map=fs5.sulc_left,axes=a,threshold=0.5)
    a.set_title('leftft > rightft, lh')

    a = ax[1][0]
    disp = plot_surf_stat_map(fs5.infl_left, dat_onsurf_lh, hemi='left', view='medial',
                              bg_map=fs5.sulc_left,axes=a,threshold=0.5)


    a = ax[0][1]
    disp = plot_surf_stat_map(fs5.infl_right,dat_onsurf_rh, hemi='right', view='lateral', 
                              bg_map=fs5.sulc_right,axes=a,threshold=0.5)
    a.set_title('leftft > rightft, rh')

    a = ax[1][1]
    disp = plot_surf_stat_map(fs5.infl_right,dat_onsurf_rh, hemi='right', view='medial',
                              bg_map=fs5.sulc_right,axes=a,threshold=0.5)



    z =  z_thrs['rightft_minus_leftft']
    dat_onsurf_lh = vol_to_surf(z, fs5.pial_left,radius=3)
    dat_onsurf_rh = vol_to_surf(z, fs5.pial_right,radius=3)

    a = ax[0][2]
    disp = plot_surf_stat_map(fs5.infl_left, dat_onsurf_lh, hemi='left', view='lateral', 
                              bg_map=fs5.sulc_left,axes=a,threshold=0.5)
    a.set_title('rightft > leftft, lh')

    a = ax[0][3]
    disp = plot_surf_stat_map(fs5.infl_left, dat_onsurf_lh, hemi='left', view='medial',
                              bg_map=fs5.sulc_left,axes=a,threshold=0.5)

    a = ax[1][2]
    disp = plot_surf_stat_map(fs5.infl_right,dat_onsurf_rh, hemi='right', view='lateral', 
                              bg_map=fs5.sulc_right,axes=a,threshold=0.5)
    a.set_title('rightft > leftft, rh')

    a = ax[1][3]
    disp = plot_surf_stat_map(fs5.infl_right,dat_onsurf_rh, hemi='right', view='medial',
                              bg_map=fs5.sulc_right,axes=a,threshold=0.5)



    z =  z_thrs['ft_minus_rest']
    dat_onsurf_lh = vol_to_surf(z, fs5.pial_left,radius=3)
    dat_onsurf_rh = vol_to_surf(z, fs5.pial_right,radius=3)

    a = ax[0][4]
    disp = plot_surf_stat_map(fs5.infl_left, dat_onsurf_lh, hemi='left', view='lateral', 
                              bg_map=fs5.sulc_left,axes=a,threshold=0.5)
    a.set_title('ft > rest, lh')

    a = ax[0][5]
    disp = plot_surf_stat_map(fs5.infl_left, dat_onsurf_lh, hemi='left', view='medial',
                              bg_map=fs5.sulc_left,axes=a,threshold=0.5)

    a = ax[1][4]
    disp = plot_surf_stat_map(fs5.infl_right,dat_onsurf_rh, hemi='right', view='lateral', 
                              bg_map=fs5.sulc_right,axes=a,threshold=0.5)
    a.set_title('ft > rest, rh')

    a = ax[1][5]
    disp = plot_surf_stat_map(fs5.infl_right,dat_onsurf_rh, hemi='right', view='medial',
                              bg_map=fs5.sulc_right,axes=a,threshold=0.5)

    outf = out_fstr + '_surfs.png'
    plt.savefig(outf)
    plt.close()


  return z_thrs,glm,evs,cntrsts,img






