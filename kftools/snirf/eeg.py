# EEG Analyses
import h5py
import mne
import numpy as np

def snirf_eeg_task_ana(f):
  
  # Load data
  
  # return a dictionary with mne object eeg data and analysis results if requested

def get_eeg(f):
    """Extracts EEG objects from a .snirf file, sets the events and channels' 10-20 montage.
    
    required input:
        fnirs2use [str]: os.path-like .snirf file address (relative or absolute)
    output:
        raw_eeg [mne.io.Raw object]: with the set events, montage, and properties
    """
    # Load H5Py file
    h5py_str = h5py.File(f, 'r')

    # Create sorted dictionary of h5py objects containing the keyword 'aux'
    foodict = sorted_nicely({k for k in h5py_str['nirs'].keys() if k.startswith('aux')})
    
    ch_names = []
    ts_eeg = np.zeros((len(foodict), np.squeeze(np.array(h5py_str['nirs']['aux16']['dataTimeSeries'])).shape[0]))
    
    # Grab the channel names & respective time series
    for xx, yy in enumerate(foodict):
      ch_names.append(np.array(h5py_str['nirs'][yy]['name']).reshape((1,1))[0][0].decode("utf-8"))
      ts_eeg[xx] = np.squeeze(np.array(h5py_str['nirs'][yy]['dataTimeSeries']))
    
    microV = [s for s in ch_names if "microV" in s]
    indices = []
    for xx in range(len(microV)):
      indices.append([i for i, x in enumerate(ch_names) if x == microV[xx]][0])
    ts_eeg_subsample = ts_eeg[indices, :]
    
    ch_names_subsample = [microV[xx][4:6] for xx in range(len(microV))]

    standard_1020_montage = mne.channels.make_standard_montage('standard_1020')
    ch_names_lower = np.array([x.lower() if isinstance(x, str) else x for x in ch_names_subsample])

    standard_1020_ch_names_lower = np.array([x.lower() if isinstance(x, str) else x for x in standard_1020_montage.ch_names]).tolist()
    for x in range(0, len(ch_names_subsample)):
        index = standard_1020_ch_names_lower.index(ch_names_lower[x])
        ch_names_subsample[x] = standard_1020_montage.ch_names[index]
    sfreq = 1000
    info = mne.create_info(ch_names=ch_names_subsample,
                           sfreq=sfreq,
                           ch_types='eeg')
    info['bads']=[]
    info['description']=f
    #Evkd data
    raw_eeg = mne.io.RawArray(ts_eeg_subsample, info)
    raw_eeg.set_montage('standard_1020')
    raw_eeg.set_eeg_reference(ref_channels='average',projection=True)
    raw_eeg.apply_proj()
    return raw_eeg  
