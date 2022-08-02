from curses import raw
import os,sys,glob,shutil,numpy as np, pandas as pd
import requests, zipfile,gdown
from datetime import datetime
from mne.io.snirf import read_raw_snirf
from kftools import data

# dictionary format : kftools -> type -> experiments -> site -> subjects
# the nifti dictionaries include 3 dictionaries each: hbo nifti, hbr nifti, and events.tsv (in this order)


def load_info(info_file=None):

  if info_file:
    f = info_file
  else:
    pth = os.path.dirname(data.__file__)
    f = os.path.join(pth,'info.txt')

  lines = [l.split(' ') for l in open(f, 'r').readlines()]
  lines = lines[1:] # ignore the column headers in the info file
  colhdrs = ['fname', 'site', 'subid', 'task', 'sesid', 'datetime', 'filetype', 'dlcode']
  rowvals = []
  for l1,l2 in lines:
    addthis = [l1] + l1.split('.')[0].split('_') + [l2.replace('\n', '')]
    rowvals.append(addthis)
  df = pd.DataFrame(rowvals, columns=colhdrs)

  return df




def fetch_file(data_dir=None,info_file=None,site='snic',task='ft',subid='sub001',sesid='ses01',
               filetype='kp-nii-evs', download_method='gdown', load_raw=True):
 
  """
  Pull selected data files and load into dictionary as mne.raw objects

  Usage:
  ------

  from kftools.data import fetch_file

  fetch_file(filetype='kp-nii-evs', sesid='ses02')

  outdir = '/external/rprshnas01/netdata_kcni/jglab/Data/kftools_data'
  fetch_file(filetype='kp-snf-hb', sesid='ses01', data_dir=outdir)


  """

  raw_dict = {} # dictionary containing mne raw objects
  
  if 'NoneType' in str(type(data_dir)): 
    data_dir = os.path.expanduser('~/.kftools')
    if not os.path.isdir(data_dir): os.makedirs(data_dir)

  df_info = load_info(info_file=info_file)
  
  qstr = ''
  for k_it,k in enumerate(['site', 'task', 'subid', 'sesid', 'filetype']):
    if 'NoneType' not in str(type(k)):
      if k_it!=0: qstr += " & "
      qstr += "%s=='%s'" %(k, eval(k))
      assert eval(k) in df_info[k].values

  idxstoget = df_info.query(qstr).index.values
  for n, idx in enumerate(idxstoget):
    dlcode,fname = df_info.loc[idx][['dlcode', 'fname']].values
    fname = os.path.join(data_dir, fname)
    if not os.path.isfile(fname):
      pull_file(dlcode,fname,download_method)
    if load_raw and 'snf' in filetype:
      raw_dict[n] = read_raw_snirf(fname)
      return raw_dict


def pull_file(dlcode,destination,download_method):

  dest_file = destination.split('/')[-1]
  print('\n\nDownloading %s\n' %dest_file)

  if download_method == "gdown":

    url = "https://drive.google.com/uc?id=" + dlcode
    gdown.download(url, destination, quiet=False)


  elif download_method == "requests":

    url = "https://docs.google.com/uc?export=download"

    session = requests.Session()
    response = session.get(
    url, params={"id": dlcode}, stream=True)

    # get the confirmation token to download large files
    token = None
    for key, value in response.cookies.items():
      if key.startswith("download_warning"):
        token = value

    if token:
      params = {"id": id, "confirm": token}
      response = session.get(url, params=params, stream=True)

    # save content to the zip-file  
    CHUNK_SIZE = 32768
    with open(destination, "wb") as f:
      for chunk in response.iter_content(CHUNK_SIZE):
        if chunk:
          f.write(chunk)




