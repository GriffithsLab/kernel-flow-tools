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

res = snirf_task_ana(f)

res.keys()


# %%
# Look at GLM results
# ---------------------------------------------------

# %% 
# HbO, Left > Right contrast

glmest = res['glm_est']
con = res['contrast_LgtR'].data
disp = kf_plot_glm_contrast_topo(glmest, con,
                                 vmin=0,vmax=5E-7,sig_thr=0.00001)

# %% 
# HbO, Right > Left contrast

glmest = res['glm_est']
con = res['contrast_RgtL'].data
disp = kf_plot_glm_contrast_topo(glmest, con,
                                 vmin=0,vmax=5E-7,sig_thr=0.00001)


# %% 
# HbR, Left > Right contrast

glmest = res['glm_est']
con = res['contrast_LgtR'].data
disp = kf_plot_glm_contrast_topo(glmest, con, chromo='hbr',
                                 vmin=0,vmax=5E-7,sig_thr=0.0001)

# %% 
# HbR, Right > Left contrast

glmest = res['glm_est']
con = res['contrast_RgtL'].data
disp = kf_plot_glm_contrast_topo(glmest, con, chromo='hbr',
                                 vmin=0,vmax=5E-7,sig_thr=0.0001)




