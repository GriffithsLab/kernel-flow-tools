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
            site='snic', task='ft', subid='sub001', sesid='ses01')

# %%
# Load the data
# ---------------------------------------------------

f = 'snic_sub001_ft_ses01_0909-1523_kp-snf-hbm.snirf'



res = load_snirf_eeg(f)

