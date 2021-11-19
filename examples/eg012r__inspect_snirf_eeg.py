# -*- coding: utf-8 -*-
r"""
=================================
Looking at KF SNIRF EEG Data
=================================


"""  

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import fetch_file
from kftools.snirf import load_snirf_eeg


# %%
# Grab the data
# --------------------------------------------------

data_dir='.'
 
fetch_file(data_dir=data_dir, filetype='kp-snf-hbm',
            site='pitch', task='ft', subid='sub010', sesid='ses01')

# %%
# Load the data
# ---------------------------------------------------

f = 'pitch_sub010_ft_ses01_1017-1706_kp-snf-hbm.snirf'
res = load_snirf_eeg(f)
res.keys()

