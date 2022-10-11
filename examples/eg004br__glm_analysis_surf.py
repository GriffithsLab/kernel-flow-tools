# -*- coding: utf-8 -*-
r"""
=================================
SNIRF HbM GLM Analysis
=================================


"""  

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import fetch_file
from kftools.snirf import snirf_task_ana,kf_plot_glm_contrast_topo


# %%
# Grab the data
# --------------------------------------------------

data_dir='.'
 
fetch_file(data_dir=data_dir, filetype='kp-snf-hbm',
            site='pitch', task='ft', subid='sub010', sesid='ses01')

# %%
# Specify and run the GLM model
# ---------------------------------------------------

f = 'pitch_sub010_ft_ses01_1017-1706_kp-snf-hbm.snirf'

res_hbo = snirf_task_ana(f,chromo='hbo')
res_hbr = snirf_task_ana(f, chromo='hbr')

res_hbo.keys()


# %%
# Look at GLM results
# ---------------------------------------------------

# %% 
# HbO, Left > Right contrast

glmest = res_hbo['glm_est']
con = res_hbo['contrast_LgtR'].data
disp = kf_plot_glm_contrast_topo(glmest, con, chromo='hbo',
                                 vmin=0,vmax=0.5,sig_thr=0.00001)
disp.show()


# %% 
# HbO, Right > Left contrast

glmest = res_hbo['glm_est']
con = res_hbo['contrast_RgtL'].data
disp = kf_plot_glm_contrast_topo(glmest, con, chromo='hbo',
                                 vmin=0,vmax=0.5,sig_thr=0.00001)
disp.show()


# %% 
# HbR, Left > Right contrast

glmest = res_hbr['glm_est']
con = res_hbr['contrast_LgtR'].data
disp = kf_plot_glm_contrast_topo(glmest, con, chromo='hbr',
                                 vmin=0,vmax=0.3,sig_thr=0.0001)
disp.show()


# %% 
# HbR, Right > Left contrast

glmest = res_hbr['glm_est']
con = res_hbr['contrast_RgtL'].data
disp = kf_plot_glm_contrast_topo(glmest, con, chromo='hbr',
                                 vmin=0,vmax=0.5,sig_thr=0.0001)
disp.show()


# %%
# On surf
from mne_nirs.visualisation import plot_glm_surface_projection
from mne.coreg import read_source_spaces,get_subjects_dir,apply_trans
from copy import copy,deepcopy
#from mne.utils import get_subjects_dir
import os
import numpy as np

subjects_dir = get_subjects_dir(raise_error=True)
fname_src_fs = os.path.join(subjects_dir, 'fsaverage', 'bem',
                                    'fsaverage-ico-5-src.fif')
src = read_source_spaces(fname_src_fs)

default_surfview_kwargs = dict(hemi = 'both',
                               size = 800,
                               figure = None,
                               clim = 'auto',
                               colormap = 'hot',
                               colorbar=True,
                               background = 'w',
                               distance=1.0,
                               mode = 'weighted',
                               surface='pial',
                               value = "test",
                               view='caudal',
                               subjects_dir = subjects_dir,
                               src = None,
                               verbose = True
                               )


glmest = res_hbo['glm_est']


con_LgtR = res_hbo['contrast_LgtR'].data
con_LgtR_t = np.nan_to_num(con_LgtR.stat())
con_LgtR_t_pos = con_LgtR_t.copy()
con_LgtR_t_pos[con_LgtR_t_pos<0] = 0

con_RgtL = res_hbo['contrast_RgtL'].data
con_RgtL_t = np.nan_to_num(con_RgtL.stat())
con_RgtL_t_pos = con_RgtL_t.copy()
con_RgtL_t_pos[con_RgtL_t_pos<0] = 0

con_RgtL_p = np.nan_to_num(con_RgtL.one_minus_pvalue())
con_LgtR_p = np.nan_to_num(con_LgtR.one_minus_pvalue())


con_LgtRest = res_hbo['contrast_LgtRest'].data
con_LgtRest_t = np.nan_to_num(con_LgtRest.stat())
con_LgtRest_t_pos = con_LgtRest_t.copy()
con_LgtRest_t_pos[con_LgtRest_t_pos<0] = 0

con_RgtRest = res_hbo['contrast_RgtRest'].data
con_RgtRest_t = np.nan_to_num(con_RgtRest.stat())
con_RgtRest_t_pos = con_RgtRest_t.copy()
con_RgtRest_t_pos[con_RgtRest_t_pos<0] = 0


inst = glmest.copy()
info = inst.copy().pick('hbo')   

df = inst.to_dataframe(order=inst.ch_names)#.iloc[idxs,:]
df = df[df.Condition=='Tapping/L']#.iloc[idxs]

k = 'LgtR_t_pos'
df[k] = con_LgtR_t_pos.copy() # con_LgtR_t_arr_nonan_nz.copy()

kws = copy(default_surfview_kwargs)
kws['value'] = k
kws['view'] = 'dorsal'#audal'
kws['mode'] = 'sum' # earest' # 'weighted'
kws['colormap'] = 'RdBu_r'
#kws['clim'] = {'kind': 'value'}#, 'lims': [-2,0.,2]}
kws['distance'] = 0.01 # 1 # 0.1 # 05
kws['src'] = src

disp = plot_glm_surface_projection(info, df, **kws)

disp.add_sensors(info.info,trans='fsaverage', 
                 fnirs=['sources', 'detectors', 'channels'])
disp.add_head(alpha=0.1)

#disp.save_image('testfig.png')
#disp.close()
#clear_output()
#Image('testfig.png')
disp.show()




