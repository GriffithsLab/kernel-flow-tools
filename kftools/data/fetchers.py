import os,sys,glob,shutil,numpy as np, pandas as pd
import requests, zipfile,gdown
from datetime import datetime
from eegnb import DATA_DIR


def fetch_dataset(
    data_dir=None,
    file_type=None,
    experiment=None,
    site="kcni",
    subjects="all",
    download_method="gdown",
    f_ids = [0]):
    """
    Return a long-form filenames list and a table saying what
    subject and session, and run each entry corresponds to
    Usage:
            data_dir = '/my_folder'
            file_type = 'snirf-mom'
            experiment = 'kf-ft'
            subjects = [1]
            # sessions = 'all'
            snirf-mom_kf-ft_fnames = fetch_dataset(data_dir=data_dir, file_type = 'snirf-mom', subjects='all', experiment='kf-ft',
            site='kcni')
            snirf-mom_rec_fnames = fetch_dataset(data_dir=data_dir, file_type = 'snirf-mom', subjects=[1], experiment='rec',
            site='kcni')
    """
    # List of file types available
    file_type_list = [
        "snirf-mom",
        "snirf-hb",
        "nifti"
    ]

    # List of experiments available
    experiments_list = [
        "kf-ft",
        "kf-fl",
        "kf-bh",
        "rec",
        "reo",
        "mmn",
    ]
  
    # List of sub available
    experiments_list = [
        "kf-ft",
        "kf-fl",
        "kf-bh",
        "rec",
        "reo",
        "mmn",
    ]

    # If no non-default top-level data path specified, use default
    if data_dir == None:
        data_dir = os.getcwd() #?

    # check parameter entries
    if experiment not in experiments_list:
        raise ValueError("experiment not in database")

    if file_type not in file_type_list:
        raise ValueError("file type not in database")

    if subjects == "all":
        subjects = [0,1,2,3,4,5,6,7,8]

    # selecting the files to download #   
    select_files = []
    if file_type == '':
      for subject in subjects:
        select_files.append(kftools_files_dict[file_type[experiment[site[subject]]]])

    # check if data has been previously downloaded
    download_it = False
    exp_dir = os.path.join(data_dir, '/kftools_datasets')
    if not os.path.isdir(exp_dir):
        download_it = True

    if download_it:
        # check if data directory exits. If not, create it
        if os.path.exists(data_dir) is not True: #% might need to change to add '/kftools_datasets'
            os.makedirs(data_dir)

        destination = os.path.join(data_dir, "/kftools_datasets") 

        if download_method == "gdown":
          for f in select_files:
            URL = "https://drive.google.com/uc?id=" + f
            gdown.download(URL, destination, quiet=False)

        elif download_method == "requests": #% can I just get rid of options for download?

            URL = "https://docs.google.com/uc?export=download"

            session = requests.Session()
            response = session.get(
                URL, params={"id": gdrive_locs[experiment]}, stream=True #%????
            )

            # get the confirmation token to download large files
            token = None
            for key, value in response.cookies.items():
                if key.startswith("download_warning"):
                    token = value

            if token:
                params = {"id": id, "confirm": token}
                response = session.get(URL, params=params, stream=True)

            # save content to the zip-file  
            CHUNK_SIZE = 32768
            with open(destination, "wb") as f:
                for chunk in response.iter_content(CHUNK_SIZE):
                    if chunk:
                        f.write(chunk)

        # unzip the file
        with zipfile.ZipFile(destination, "r") as zip_ref:
            zip_ref.extractall(data_dir)

        # # remove the compressed zip archive
        # os.remove(destination)

    if subjects == "all":
        subjects = [0,1,2,3,4,5,6,7,8]

    # if sessions == "all":
    #     sessions = ["*"]

    # If 'all' subjects and 'all' sessions:
    # if (subjects[0] == "*") and (sessions[0] == "*"):
    #     pth = os.path.join(
    #         exp_dir, f"subject{subjects[0]}", f"session{sessions[0]}", "*.csv"
    #     )
    #     fnames = glob.glob(pth)

    # Else, if specific subjects and sessions
    else:
        fnames = []
        for subject_nb in subjects:
            if subject_nb != "*":
                # Format to get 4 digit number, e.g. 0004
                subject_nb = float(subject_nb)
                subject_nb = "%03.f" % subject_nb
                # for session_nb in sessions:
                #     # Formt to get 3 digit number, e.g. 003
                #     if session_nb != "*":
                #         session_nb = float(session_nb)
                #         session_nb = "%02.f" % session_nb
                pth = os.path.join(exp_dir,f"subject{subject_nb}",
                  f"session{session_nb}","*.csv" )
        fpaths = glob.glob(pth)
        fnames += fpaths

    return fnames



# def zip_data_folders(experiment: str,
#                      site: str="local_ntcs"):

#     """
#     Run data zipping
#     Usage
#     from eegnb.datasets.datasets import zip_data_folders
#     zip_data_folders(experiment='visual-N170')
#     See also the command-line program
#     eegnb runzip -ip
#     """

#     print('\nRunning Data Zipper')
#     zip_directory=os.path.join(DATA_DIR,experiment,site)
#     print('Looking for {} within {} \n'.format(experiment+'/'+site,DATA_DIR))
    
#     if not os.path.isdir(zip_directory):
#         print('Invalid Directory')
#         raise ValueError ('{} directory does not exist'.format(zip_directory))

#     print('Files Found! Zipping all files in {} '.format(zip_directory))

#     date_time=datetime.now()
#     datetime_str=date_time.strftime("%d_%m_%Y_%H:%M")
#     output_filename=os.path.join(os.path.expanduser("~/Desktop"),experiment+'_'+site+'-'+datetime_str+'_zipped')
    
#     shutil.make_archive(output_filename,'zip',zip_directory)
#     print('Zip file location is at {}\n '.format(output_filename))



# dictionary format : kftools -> type -> experiments -> site -> subjects
# the nifti dictionaries include 3 ids each: hbo nifti, hbr nifti, and events.tsv (in this order)
kftools_files_dict = {
    'snirf-mom': {
        'kf-ft':{
            'kcni':{
                0:['11Lcv-dczig_SLGHcyl8hstE3IS5R8oXO'],
                1:[],
                2:[],
                3:[],
                4:[],
                5:[],
                6:[],
                7:[],
                8:[]
            },
            'advc':{
                0:[],
                1:[],
                2:[],
                3:[],
                4:[],
                5:[],
                6:[],
                7:[],
                8:[]
            },
            'kern':{
                0:[],
                1:[],
                2:[],
                3:[],
                4:[],
                5:[],
                6:[],
                7:[],
                8:[]
            }
        },
        'kf-fl':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'kf-bh':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'rec':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'reo':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'mmn':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
    'snirf-hb': {
        'kf-ft':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7:[18Zkh-NQtPRHv_oWBFjB9QUg5jMGPidlK],
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'kf-fl':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'kf-bh':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'rec':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'reo':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'mmn':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        }
    'nifti': {
        'kf-ft':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7:[15uLx8vGhLWyY_Wib-X3Bd3VAzd-I7vU9,15vwsTYcSOVJf2qkmGt8V6yMQ8oLcGJDz,18tY67pVCQbdQPv_AGOCuxK4Rsp3QrbFY)],
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'kf-fl':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'kf-bh':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'rec':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'reo':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        },
        'mmn':{
            'kcni':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'advc':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            },
            'kern':{
                0:,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8
            }
        }

}

