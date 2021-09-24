# KernelFlow Dev Sprint 

## Overview / Goals: What is Kernel Flow, why is it important, and what are we doing? 
* In the area of non-invasive optical brain imaging devices and techniques, time-domain near-infrared spectroscopy (TD-NIRS) stands as a more effective method for capturing and analyzing data (when compared to continuous wave (CW) systems for its ability to be analyzed with a greater content of information). With the lack of commercial systems available, adoption of the technology has been slow (Ban, et al., 2021). 
* Using this technology, we can collect functional near-infrared spectroscopy (fNIRS) data (with much higher resolution than functional magnetic resonance imaging (fMRI)). This type of data and work has become a reference tool and emerging paradigm in research on experimental motor control (Guerin et al., 2021), shown use in studying functional connectivity (FC) structures of the brain in resting states and at sleep states (Nguyen et al., 2018), ...
* During this sprint, we’ve attempted to create a procedure of experimentation (obtaining and analyzing results of individuals performing specific tasks using the Kernel Flow software) and theoretical analysis (processing, visualizing, and aggregating/summarizing results of the experiments) for the Kernel Flow technology with the hopes of establishing its use among labs. We’re focusing specifically on TD-fNIRS (functional near-infrared spectroscopy) 
* During this sprint, we’re working with researchers across the world together to create and establish a procedure for collecting data and analyzing it through the Flow TD-FNIRS. 
 
## Kernel Flow experimentation: Taha Morshedzadeh, Shreyas Harita, Parsa Oveisi, Renee Huang
### Protocol / procedure:
* Setup helmet 
* Take head measurements
* Example measurements:
1. Nasion to Inion - 40 cm (for Shrey), 36 cm (for Davide)
2. Ear to Ear - 40 cm (for Shrey), 38 cm (for Davide)
3. Head Circumference - 56 cm (for Shrey), 59 cm (for Davide).
* Begin Kernel Flow application: 
* Initialize the experiment
* Tune lasers
* View the EEG coupling

### Record / collect the data
* For each of the following, collect fNIRS and EEG data in a time-synchronized fashion. Use kernel acquisition method, 
1. Dummy source analysis
2. Finger tapping
3. Visualizing checkerboard data (already completed from the day before)
4. Breath-holding tasks

## Python-based analysis (mne-nirs): Davide Momi, Frank Mazza, FuTe Wong, Andrew Clappison, Kevin Kadak, John Griffiths, Jerry Jeyachandra
* Using MNE for visualizing the data obtained from the experiments
* Pre-processing and data analysis
* Event timing reproduced from the tasks 
* Script created to fix the irregular sampling in a SNIRF file (jg_fix_irregular_sampling_in_snirf_file.ipynb) 
* Topoplots produced 

## Python-based analysis (Nilearn / Nifti files): Shreyas Harita, Kevin Kadak, Hussain Ather
* Extract and visualize the `nifti` file data using Nilearn

## Matlab-based analysis (homer3): Sorenza, Zheng, Hussain, Colleen, Andreea, Gabrielle
* Used Homer3 (and AtlasViewer) for analyzing fNIRS data to obtain estimates and maps of brain activation that can be visualized with AtlasViewer (for mapping the data onto various parts of the brain). 
* Run through steps for Homer3 basic analysis with example data
* Visualize data with AtlasViewer 
* Transform the 2-D probe into a 3-D format
* This involved using the information specific to the probe based on how the 2-D layout corresponds to the 3-D layout of the brain itself. Use the source distance and relative positions of the optodes. 
* We also used the anchor points, dummy optodes, and springs in Homer3 and the AtlasViewer GUI. 
* Load the kernel .snirf file into Homer3.
* Found that the .snirf file created by the kernel does not mention the (necessary) wavelength index.
* AtlasViewer also offers ways of letting one locate the probe on the scalp so that we can actually measure our target brain region beneath, project and visualize HbO/HbR responses onto the brain surface, obtain the MNI coordinates of the optode/channel locations, and test the probe localization and fabrication error AtlasViewer. 
* Used the scripts in the Homer3 MATLAB toolbox for extracting and displaying information of the .snirf files themselves. 
* The `kernelsprintday1.m` script gives one way of extracting and visualizing the data. 
* Trying to output volume data to nifti files
* Using NIfTI and ANALYZE tools

### References

Ban, Han Y., et al. "Kernel flow: a high channel count scalable TD-fNIRS system." Integrated Sensors for Biological and Neural Sensing. Vol. 11663. International Society for Optics and Photonics, 2021. 

Guérin, Ségolène MR, et al. "Effects of Motor Tempo on Frontal Brain Activity: An fNIRS Study." NeuroImage 230 (2021): 117597.

Nguyen, Thien, et al. "Exploring brain functional connectivity in rest and sleep states: a fNIRS study." Scientific reports 8.1 (2018): 1-10.
