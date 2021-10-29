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


# %%
# Grab the data
# --------------------------------------------------

f1a = 'pitch_sub010_ft_ses01_1017-1706_kp-nii-hbo.nii.gz'

f1b = 'pitch_sub010_ft_ses01_1017-1706_kp-nii-hbr.nii.gz'

f1c = 'pitch_sub010_ft_ses01_1017-1706_kp-nii-evs.tsv'
    
res =  ft_glm_ana(f1a, f1b, out_fstr = 'test_')#
z_thrs,glm,evs,cntrsts,img = res




