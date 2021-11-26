# Homer README
Malab codes

# Instructions to run Homer3 and AtlasViewer pipeline script in SCC:

1)	Go to folder location with probe.SD, digpts.txt, Homer3, AtlasViewer and snirf  subject data file
2)	Load MATLAB modules
  -	module load math/MATLAB/2020b
  -	Open MATLAB: matlab –nodesktop
3)	Steps before running pipeline script
  -	Add Homer3 to path: cd Homer3  setpaths
  -	Add AtlasViewer to path: cd ..  avpath = fullfile(pwd, ‘AtlasViewer’))  addpath(genpath(avpath)) 
  -	Set name of AtlasViewer path: dirname_atlas_viewer = fullfile(pwd,avpath)
  -	Set current folder path: currFolder = pwd
  -	Set snirf file name: snirfFileName = fullfile(pwd, ‘insert_name_file.snirf’)
4)	Run pipeline code
  -	[dc, fwmodel] = Homer_Atlas(snirfFileName, dirname_atlas_viewer, currFolder)
 
# What does the pipeline code do?

1)	[fwmodel, probe] = ForwardModel(dirname_atlas_viewer, currFolder)
  -	Run ForwardModel function from AtlasViewer to get the sensitivity profile
2)	MeanMoments(snirfFileName): From the original snirf data subject file, the function extracts the mean moments
3)	HomerProcessing(snirf_data): 
  -	hmrR_PreprocessIntensity_NAN(res.data): Remove all NAN values and replace with the mean moments
  -	hmrR_Intensity2OD(acquired.data): Intensity to change in optical density
  -	hmR_BandpassFilt(dod_nan_rm, 0.01, 0.50): Apply bandpass filter to data
  -	hmrR_OD2Conc_zw(dod, acquired.probe, probe, [1,1]): Converts optical density to concentration for existing modules
  -	SaveData(dc, fwmodel.Adot): Function that saves HbO and HbR time series and their corresponding channel to vertice transforms used in python code

