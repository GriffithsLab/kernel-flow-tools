

from kftools.data import fetch_file

data_dir = '/external/rprshnas01/netdata_kcni/jglab/Data/kftools_data'


# Nii Events file
fetch_file(data_dir=data_dir, filetype='kp-nii-evs')

# Nii HbO file
fetch_file(data_dir=data_dir, filetype='kp-nii-hbo')

# HB Moments file
fetch_file(data_dir=data_dir, filetype='kp-snf-hbm')



