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
# HbO analyses
# --------------------------------------------------    

# %%
# Fit model
res_hbo =  ft_glm_ana(nii_hbo_f, nii_ev_f, out_fstr = '')
z_thrs,glm,evs,cntrsts,img = res_hbo

# %%
# Glass brain plots
z = z_thrs['rightft_minus_leftft']
disp = plot_glass_brain(z,colorbar=True, threshold=1)#,vmin=10,vmax=20);
a.set_title('right > left');

z = z_thrs['leftft_minus_rightft']
disp = plot_glass_brain(z,colorbar=True, threshold=1)#,vmin=10,vmax=20);
a.set_title('left > right');



# %% 
# HbR analyses
# -------------------------------------------------

res_hbr =  ft_glm_ana(nii_hbr_f, nii_ev_f, out_fstr = '')
z_thrs,glm,evs,cntrsts,img = res_hbr

# %%
# Glass brain plots
z = z_thrs['rightft_minus_leftft']
disp = plot_glass_brain(z,colorbar=True, threshold=1)#,vmin=10,vmax=20);
a.set_title('right > left');

z = z_thrs['leftft_minus_rightft']
disp = plot_glass_brain(z,colorbar=True, threshold=1)#,vmin=10,vmax=20);
a.set_title('left > right');



