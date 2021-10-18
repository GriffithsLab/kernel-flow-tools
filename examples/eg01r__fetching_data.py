# -*- coding: utf-8 -*-
r"""
=================================
Downloading KFTools Example Data
=================================


"""  # noqa:E501

from kftools.data import fetch_file


#data_dir = '/external/rprshnas01/netdata_kcni/jglab/Data/kftools_data'
data_dir=None

# sphinx_gallery_thumbnail_number = 1

# %%
# First we do a thing etc etc.


# %%
# Nii Events File
# ---------------------------------------------------
#
fetch_file(data_dir=data_dir, filetype='kp-nii-evs')

# %%
# Nii HbO file
# ---------------------------------------------------
fetch_file(data_dir=data_dir, filetype='kp-nii-hbo')

# %%
# HB Moments file
# ---------------------------------------------------
fetch_file(data_dir=data_dir, filetype='kp-snf-hbm')




