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


##### create brain mask for sub-control035

cd $new_fs/sub-control035/SUMA

# rename fs_parc_wb_mask.nii.gz
3dcopy fs_parc_wb_mask.nii.gz fs_parc_wb_mask_orig.nii.gz
rm fs_parc_wb_mask.nii.gz

# copied from adjunct_suma_fs_mask_and_qc

# make a filled mask from aparc+aseg (output to ../. directory)
3dmask_tool                                                       \
    -overwrite                                                    \
    -dilate_inputs 2                                              \
    -prefix mask1.nii.gz                                          \
    -inputs aseg.auto.nii.* 

3dmask_tool                                                       \
    -overwrite                                                    \
    -fill_holes                                                   \
    -prefix mask2.nii.gz                                          \
    -inputs mask1.nii.gz

3dmask_tool                                                       \
    -overwrite                                                    \
    -datum         byte                                           \
    -dilate_inputs -2                                             \
    -prefix mask3.nii.gz                                          \
    -inputs mask2.nii.gz

3dcopy                                                            \
    -overwrite                                                    \
    mask3.nii.gz                                                  \
    fs_parc_wb_mask.nii.gz

