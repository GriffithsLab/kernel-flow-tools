# -*- coding: utf-8 -*-
r"""
=================================
SNIRF EEG MNE Task Analysis
=================================


"""  

# sphinx_gallery_thumbnail_number = 1

# %%
# Importage
# --------------------------------------------------

# KF Tools and related imports
from kftools.data import fetch_file
#from kftools.snirf import snirf_eeg_task_ana


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

#res = snirf_eeg_task_ana(f)

#res.keys()


