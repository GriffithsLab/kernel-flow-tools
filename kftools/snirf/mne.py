

import os,sys,glob,numpy as np,pandas as pd

from matplotlib import pyplot as plt

import mne
from mne.io import read_raw_snirf

from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm
import mne_nirs

from nilearn.plotting import plot_design_matrix
from mne_nirs.statistics import run_glm

from mne.viz.utils import _plot_sensors,_check_sphere
from mne.viz import plot_sensors

from mne.viz.topomap import plot_topomap
from mne_nirs.visualisation._plot_GLM_topo import _plot_glm_topo


import matplotlib as mpl
import mne
from copy import deepcopy
from mne_nirs.visualisation._plot_GLM_topo import _handle_overlaps



import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

import mne
from mne import Info
from mne.utils import warn
from mne.channels.layout import _merge_ch_data
from mne.io.pick import _picks_to_idx, _get_channel_types





def run_pipeline(f,subselect_with = None, subselect_range = None):

  hbm = read_raw_snirf(f, preload=True)#, optode_frame='unknown', preload=False, verbose=None)#[source]¶

  hbm.resample(sfreq=1) 

  #hbm._data = np.nan_to_num(hbm._data)

  dropchans = [c for c in hbm.ch_names if 'HbR' in c]
  hbo = hbm.copy()
  hbo.drop_channels(dropchans);

  hbo.set_channel_types({chn: 'hbo' for chn in hbo.ch_names })

  for c_it,c in enumerate(hbo.ch_names): 
    hbo.info['chs'][c_it]['loc'][:3] = hbo.info['chs'][c_it]['loc'][6:9]


  evs = get_events_from_snirf(f)
  for e_it, e in enumerate(evs['BlockType']):
    if ((e == 'Left') or (e =='Right')):
      txt = hbo.annotations.description[e_it] 
      txt = txt.replace('StartTrial', 'Tapping') + '/' + e[0]
      hbo.annotations.description[e_it] = txt
    dur = evs.iloc[e_it].loc['Duration']
    hbo.annotations.duration[e_it] = dur

  events, event_dict = mne.events_from_annotations(hbo)

  """
  equal_idx = []
  for xx in range(events.shape[0]-1):
    if events[xx,0] == events[xx+1,0]:
      equal_idx.append(xx)
  equal_idx.reverse()
  events = np.delete(events, equal_idx, 0)
  """;

  event_dict.pop('StartExperiment', None)
  event_dict.pop('StartBlock/L', None)
  event_dict.pop('StartBlock/R', None)
  event_dict.pop('StartIti/L', None)
  event_dict.pop('StartIti/R', None)

  hbo.annotations.delete(hbo.annotations.description=='StartExperiment')
  hbo.annotations.delete(hbo.annotations.description=='StartBlock/L')
  hbo.annotations.delete(hbo.annotations.description=='StartBlock/R')
  hbo.annotations.delete(hbo.annotations.description=='StartIti/L')
  hbo.annotations.delete(hbo.annotations.description=='StartIti/R')


  tmin = -5 # -0.1
  tmax = 20 # 0.5
  hbo_ep = mne.Epochs(hbo, events, event_id=event_dict,
                      tmin=tmin, tmax=tmax,
                      reject_by_annotation=True,
                      proj=False,#True, 
                      baseline=(0, 0), 
                      preload=True,
                      detrend=None, verbose=True,
                      event_repeated='merge')

  hbo_ep._data = np.nan_to_num(hbo_ep._data)

  hbo_ev_L = hbo_ep['Tapping/L'].average()
  hbo_ev_R = hbo_ep['Tapping/R'].average()

  hbo_ev_LgtR = mne.combine_evoked([hbo_ev_L, hbo_ev_R], weights=[1, -1])
  hbo_ev_RgtL = mne.combine_evoked([hbo_ev_R, hbo_ev_L], weights=[1, -1])


  s = mne_nirs.experimental_design.create_boxcar(hbo)
  design_matrix = make_first_level_design_matrix(hbo,
                                               drift_model='cosine',
                                               high_pass=0.005,  # Must be specified per experiment
                                               hrf_model='spm',
                                               stim_dur=5.0)
  # Fit to sphere
  # (dug this out from the plot_sensors function calls)
  sphere_params =  _check_sphere('auto', info=hbo.info, sphere_units='m')

  if subselect_with == 'notnan':
      pickthese = np.nonzero(np.isnan(hbo._data.sum(axis=1))==False)[0]
  elif type(subselect_range) == list:
      pickthese = range(subselect_range[0], subselect_range[1])
  else:
      pickthese = hbo.ch_names    
  data_subset = hbo.copy().pick(picks=pickthese)
  glm_est = run_glm(data_subset, design_matrix,n_jobs=5)

  contrast_matrix = np.eye(design_matrix.shape[1])
  basic_conts = dict([(column, contrast_matrix[i])
                    for i, column in enumerate(design_matrix.columns)])
  
  contrast_vec_LgtR = basic_conts['Tapping/L'] - basic_conts['Tapping/R']
  contrast_LgtR = glm_est.compute_contrast(contrast_vec_LgtR, contrast_type='t')
  #contrast_LgtR.plot_topo(sphere='auto');#,vmin=-1E4,vmax=1E4);
  estimates = contrast_LgtR.data.effect[0]
  info = contrast_LgtR.info
  estmrg_LgtR, pos_LgtR, chs_LgtR, sphere_LgtR = _handle_overlaps(info, t,
                                                                  sphere_params, estimates)

 
  contrast_vec_RgtL = basic_conts['Tapping/R'] - basic_conts['Tapping/L']
  contrast_RgtL = glm_est.compute_contrast(contrast_vec_RgtL, contrast_type='t')
  #contrast_RgtL.plot_topo(sphere='auto')#,vmin=-1E4,vmax=1E4);
  estimates = contrast_RgtL.data.effect[0]
  info = contrast_RgtL.info
  estmrg_RgtL, pos_RgtL, chs_RgtL, sphere_RgtL = _handle_overlaps(info, t,
                                                                  sphere_params, estimates)

  returnstuff = dict(evs=event_dict,hbo=hbo,hbo_ep=hbo_ep,
                     hbo_ev_LgtR=hbo_ev_LgtR,
                     hbo_ev_RgtL=hbo_ev_RgtL,
                     s=s,design_matrix=design_matrix,
                     data_subset=data_subset,
                     sphere_params=sphere_params,
                     glm_est=glm_est,
                     contrast_LgtR = contrast_LgtR,
                     contrast_vec_LgtR = contrast_vec_LgtR,
                     contrast_RgtL = contrast_RgtL,
                     estmrg_LgtR = estmrg_LgtR, 
                     estmrg_RgtL = estmrg_RgtL,
                     pos_LgtR = pos_LgtR,
                     pos_RgtL = pos_RgtL, 
                     chs_LgtR = chs_LgtR,
                     chs_RgtL = chs_RgtL,
                     sphere_LgtR = sphere_LgtR,
                     sphere_RgtL = sphere_RgtL)


  return returnstuff






def plot_glm_contrast_topo(inst, contrast, figsize=(12, 7), sphere=None,
                              vmin=None,vmax=None,positive_only=True,cmap=False):

    info = deepcopy(inst if isinstance(inst, Info) else inst.info)

    # Extract types. One subplot is created per type (hbo/hbr)
    types = np.unique(_get_channel_types(info))

    # Extract values to plot and rescale to uM
    estimates = contrast.effect    #[0]
    
    # estimates = estimates * 1e6

    estimates = np.squeeze(estimates) # JG_ADD

    if positive_only: estimates[estimates<0] = 0 # JG_ADD

    # Create subplots for figures
    fig, axes = plt.subplots(nrows=1,
                             ncols=len(types),
                             figsize=figsize)
    # Create limits for colorbar
    if not vmax:
      vmax = np.max(np.abs(estimates))
    if not vmin:
      vmin = vmax * -1.
    if not cmap:
      cmap = mpl.cm.RdBu_r
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    for t_idx, t in enumerate(types):


        estmrg, pos, chs, sphere = _handle_overlaps(info, t, sphere, estimates)


        # Deal with case when only a single chroma is available
        if len(types) == 1:
            ax = axes
        else:
            ax = axes[t_idx]

        # Plot the topomap
        mne.viz.topomap.plot_topomap(estmrg, pos,
                                     extrapolate='local',
                                     names=chs,
                                     vmin=vmin,
                                     vmax=vmax,
                                     cmap=cmap,
                                     axes=ax,
                                     show=False,
                                     sphere=sphere)
        # Sets axes title
        if t == 'hbo':
            ax.set_title('Oxyhaemoglobin')
        elif t == 'hbr':
            ax.set_title('Deoxyhaemoglobin')
        else:
            ax.set_title(t)

    # Create a single colorbar for all types based on limits above
    ax1_divider = make_axes_locatable(ax)
    cax1 = ax1_divider.append_axes("right", size="7%", pad="2%")
    cbar = mpl.colorbar.ColorbarBase(cax1, cmap=cmap, norm=norm,
                                     orientation='vertical')
    cbar.set_label('Contrast Effect', rotation=270)

    return fig



