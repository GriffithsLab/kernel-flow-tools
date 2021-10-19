# -*- coding: utf-8 -*-
r"""
=================================
Downloading KFTools Example Data
=================================


"""  

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import load_info,fetch_file
from mne.io import read_raw_snirf

# Some neuroimaging library things
import nibabel as nib
from nilearn.plotting import plot_stat_map
from nilearn.image import index_img

# Some generic things
from matplotlib import pyplot as plt
import os
import pandas as pd


# %%
# Specify data download location
# --------------------------------------------------

# Give a specific location 
# data_dir = '/external/rprshnas01/netdata_kcni/jglab/Data/kftools_data'
# or maybe
# data_dir = '.'

# If data_dir=None the default location is used, which is '~/.kftools' 
data_dir=None


# %%
# List available files
# --------------------------------------------------
info = load_info()
info[['fname', 'site','subid', 'task', 'sesid', 'datetime', 'filetype']]


# %%
# HB Moments file
# ---------------------------------------------------

# The kernel portal gives two .snirf file options: 
# 1. The 'raw' moments file
# 2. A 'Hb moments' file

# The Hb moments file has some initial preprocessing applied to it, 
# including optical density and modified beer lambert law calculations.
 
fetch_file(data_dir=data_dir, filetype='kp-snf-hbm',
           site='snic', task='ft', sub='sub001', ses='ses01')

raw = read_raw_snirf('~/.kftools/snic_sub001_ft_ses01_0909-1523_kp-snf-hbm.snirf')

df_raw = raw.to_data_frame()

fig, ax = plt.subplots(figsize=(12,3))
df_raw[raw.ch_names[0:5]].loc[3000:].plot(ax=ax)


# %%
# Nifti file
# ---------------------------------------------------
#

# In-the-brain, DOT-reconstructed HbO and HbR time series 
# are provided in the form of 4D nifti images. 

fetch_file(data_dir=data_dir, filetype='kp-nii-hbo')
           site='snic', task='ft', sub='sub001', ses='ses01')

f = '~/.kftools/snic_sub001_ft_ses01_0909-1523_kp-nii-hbo.nii.gz'
imgs = nib.load(f)
img = index_img(imgs,0)
disp = plot_stat_map(img)


# Note that for task-based statistical models, there is a 
# separate and modified events file, that contains the same 
# information as can be pulled directly from the .snirf files
# (see later examples), but with modified timings. 

fetch_file(data_dir=data_dir, filetype='kp-nii-evs',
           site='snic', task='ft', sub='sub001', ses='ses01')

f = '~/.kftools/snic_sub001_ft_ses01_0909-1523_kp-nii-evs.tsv'
df = pd.read_csv(f, sep='\t')
df


