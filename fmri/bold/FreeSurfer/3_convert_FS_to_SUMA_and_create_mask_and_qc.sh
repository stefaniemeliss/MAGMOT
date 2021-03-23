#!/bin/bash

################################################################################
# converting FreeSurfer to AFNI to create masks for pre-processing
################################################################################


# this script has to be run after the FS recon-all command
module load freesurfer6.0.0
module load afni19.3.03

# define DIR
DIR=/storage/shared/research/cinn/2018/MAGMOT
BIDS_dir=$DIR/rawdata
deriv_dir=$DIR/derivatives
old_fs=$deriv_dir/FreeSurfer
new_fs=$deriv_dir/freesurfer
mkdir $new_fs

# change directory to the raw NIFI files
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
#subjects=(sub-experimental005)
#subjects=(sub-control003)

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

	##### copy FreeSurfer output to new directory ##### 

	#determine folders
	old_subj=$old_fs/$subject
	new_subj=$new_fs/$subject

	# make directories
	mkdir $new_subj

	# copy folders
	cp -r $old_subj/surf  $new_subj/surf
	cp -r $old_subj/mri  $new_subj/mri

	# cd into subject folder
	cd $new_subj

	##### use SUMA to convert FreeSurfer output ##### 
	echo START WITH SUMA CONVERSION
	@SUMA_Make_Spec_FS -sid ${subject} -NIFTI

	# define SUMA folder
	suma_dir=$new_subj/SUMA

	##### create mask and QC images #####
	echo CREATE THE WM AND VENTRICLE MASK 
	adjunct_suma_fs_mask_and_qc -sid ${subject} -suma_dir $suma_dir

done

