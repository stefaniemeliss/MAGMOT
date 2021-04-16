#!/usr/bin/python

# The osfclient is a python library and a command-line client for up- and downloading files to and from your Open Science Framework projects.

# load module
module load anaconda3

# install osfclient
#pip install --upgrade pip
#pip install six==1.14.0
#pip install osfclient

# determine OSF project id
projectid="fhqb7"

# upload a single file to an OSF project
#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/pyfMRIqc/pyfMRIqc_output.csv osfstorage/data/fmri/pyfMRIqc_output.csv

# upload data from afni quality check to OSF
#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/extents.txt osfstorage/data/fmri/extents.txt

#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/motion.txt osfstorage/data/fmri/motion.txt
#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/TSNR.txt osfstorage/data/fmri/TSNR.txt

# upload files to MMC data directory
projectid="eyzwb"

#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/pyfMRIqc/pyfMRIqc_output.csv osfstorage/data_pyfMRIqc.csv

# upload data from afni quality check to OSF
#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/data_extents.txt osfstorage/data_extents.txt
#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/data_motion.txt osfstorage/data_motion.txt
#osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/data_TSNR.txt osfstorage/data_TSNR.txt
osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force /storage/shared/research/cinn/2018/MAGMOT/derivatives/analysis/MMC_paper/sFC/data_sFC.txt osfstorage/data_sFC.txt

