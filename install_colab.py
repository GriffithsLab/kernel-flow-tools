import os

# Replace py38 with py37 in setup file 
txt = open('setup.py', 'r').read()
open('setup_orig.py', 'w+').write(txt)
txt_py37 = txt.replace('3.8', '3.7')
open('setup.py', 'w').write(txt_py37)

# Install repo
os.system('pip intall -r requirements.txt')
os.system('pip install -e .')
