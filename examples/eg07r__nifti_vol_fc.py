# -*- coding: utf-8 -*-
r"""
=================================
Nifti volumes - functional connectivity analyses
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

"""

### seed coords from (Eggebrecht et al., 2014) ###

seed_coord_dict = {'vis':(-19.5, -102, -3), 'aud':(-67.5, -27, 12), 'mot':(-67.5, -12, 27), 'DAN':(-58.5, -69, -6), 'FPC':(-52.5, 24, 33), 'DMN':(-43.5, 21, 51)}
seed_coord_ser = pd.Series(seed_coord_dict)

### making MNI mask ###

# import nilearn
# nilearn.__file__

avg152_f = '/usr/local/lib/python3.7/dist-packages/nilearn/datasets/data/avg152T1_brain.nii.gz'
avg152_img = nib.load(avg152_f)
avg152_dat = avg152_img.get_data()
avg152_mask_dat = (avg152_dat > 0).astype(float)
avg152_mask_img = nib.Nifti1Image(avg152_mask_dat, affine=avg152_img.affine)


nii_fs = glob.glob(fromkcni_data_dir + '/*rec*hbo*.nii*')
f = nii_fs[0]
k = 'DMN'
z_maps,z_maps_masked = fc_for_seeds(f,dothese=[k])#,clip_vols=[10,10])
img = z_maps[k]
img_thr,_ = threshold_stats_img(img,threshold=thr,cluster_threshold=50,alpha=0.05)
disp = plot_glass_brain(img_thr,plot_abs=False,colorbar=True);
disp.add_markers(marker_coords=[seed_coord_dict[k]], marker_color='g', marker_size=300)

""";
