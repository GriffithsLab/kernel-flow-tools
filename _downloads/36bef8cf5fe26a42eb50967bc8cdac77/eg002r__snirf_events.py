# -*- coding: utf-8 -*-
r"""
=================================
Getting events from .snirf files
=================================


"""  

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import fetch_file
from kftools.snirf import get_events_from_snirf
from mne.io import read_raw_snirf


# %%
# Grab the data
# --------------------------------------------------

data_dir='.'
 
fetch_file(data_dir=data_dir, filetype='kp-snf-hbm',
           site='snic', task='ft', subid='sub001', sesid='ses02')


# %%
# Get the events
# --------------------------------------------------

f = 'snic_sub001_ft_ses02_0917-1313_kp-snf-hbm.snirf'
evs = get_events_from_snirf(f)
evs


