#!/usr/bin/bash

# load anaconda module
module load anaconda3

# change working directory to where the script is saved
cd /storage/shared/research/cinn/2018/MAGMOT/scripts/fmri/bold/cut/

# call python script
python MAGMOT_cut_before_preprocessing.py
