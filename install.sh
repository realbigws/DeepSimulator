#!/bin/bash


#-> 1. create an conda environment
conda remove --name tensorflow_cdpm --all -y
conda create --name tensorflow_cdpm python=2.7 -y

#-> 2. activate the environment, install the following required modules
source activate tensorflow_cdpm
pip install tensorflow==1.2.1
pip install tflearn==0.3.2
pip install tqdm==4.19.4
pip install scipy==0.18.1
pip install h5py==2.7.1
pip install numpy==1.13.1
pip install sklearn
source deactivate

