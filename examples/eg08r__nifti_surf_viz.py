# -*- coding: utf-8 -*-
r"""
=================================
Nifti-derived surfaces - visualization
=================================

Blah blah

blah
blah
blah


*****************
Thing
*****************

Thingywot
==============

Blah blah


***********************
Anotherthing
***********************


Whatevs
=======================

Stuff
Stuff
Stuff

"""  # noqa:E501

import numpy as np
import pandas as pd
import mne

# sphinx_gallery_thumbnail_number = 1

# %%
# First we do a thing etc etc.

pd.DataFrame(np.random.normal(size=(16, 100))).plot(kind='hist')

# %%
#
# .. warning:: athing a not etc
#              blah blah
thisvar  = ['hbo', 'hbr', 'hbo', 'hbr',
            'hbo', 'hbr', 'hbo', 'hbr']

# %%
# More stuff etc blah
# blah blah
montage = mne.channels.make_standard_montage('artinis-octamon')
#raw.set_montage(montage)

# View the position of optodes in 2D to confirm the positions are correct.
print('thing')#           ``something here`` to ``'eg'``.

#
# More stuff
# More stuff
#
# .. note:: Things


# %%
# Applying standard sensor locations to imported data
# ---------------------------------------------------
#
# Having information about optode locations may assist in your analysis.
# Beyond the general benefits this provides (e.g. creating regions of interest,
# etc), this is may be particularly important for fNIRS as information about
# the optode locations is required to convert the optical density data in to an
# estimate of the haemoglobin concentrations.
# MNE-Python provides methods to load standard sensor configurations
# (montages) from some vendors, and this is demonstrated below.
# Some handy tutorials for understanding sensor locations, coordinate systems,
# and how to store and view this information in MNE-Python are:
# :ref:`tut-sensor-locations`, :ref:`plot_source_alignment`, and
# :ref:`ex-eeg-on-scalp`.
#
# Below is an example of how to load the optode positions for an Artinis
# OctaMon device.
#
# .. note:: It is also possible to create a custom montage from a file for
#           fNIRS with :func:`mne.channels.read_custom_montage` by setting
#           ``coord_frame`` to ``'mri'``.

print(5)

# View the position of optodes in 2D to confirm the positions are correct.
print(6)
