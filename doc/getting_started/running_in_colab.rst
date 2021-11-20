======================================
Running kftools in google colab
======================================

Google colab is an increasingly convenient and useful scientific omputing environment. 

KFTools can be used in google colab with some minimal additional steps. 

First, in a browser tab, navigate to https://colab.research.google.com/github/GriffithsLab/kernel-flow-tools , 
selecting organization GriffithsLab, repository kernel-flow-tools, and branch gh-pages. 

You should then see a drop-down list of .ipynb files, corresponding to each of the examples in the sphinx-gallery docs. 

Select one of these, and approve the 'run anyway' option. 

Then, insert and run the following in a cell at the top:


Clone latest version from github

.. code::

    $ import os,time
    $ os.system('rm -rf kernel-flow-tools')
    $ os.system('git clone https://github.com/griffithslab/kernel-flow-tools')
    $ time.sleep(3)
    $ os.chdir('kernel-flow-tools')
    $ time.sleep(3)
    $ os.system('python install_colab.py')    
    
Now you should be good to continue with the rest of the example code in the notebook, and experiment with new ideas. 

