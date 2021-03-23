#!/bin/bash

################################################################################
# converting FreeSurfer to AFNI to create masks for pre-processing
################################################################################


# this script has to be run after the FS recon-all command
module load freesurfer6.0.0
module load afni19.3.03


# define DIR
DIR=/storage/shared/research/cinn/2018/MAGMOT

# change directory to the raw NIFI files
cd $DIR/MAGMOT_BIDS/

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
#subjects=(sub-experimental005)
#subjects=(sub-control003)

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

	# define BIDS anat folder
	anat_dir=$DIR/MAGMOT_BIDS/"${subject}"/anat
	anat_dir=$DIR/derivatives/afniproc/"${subject}"/SSwarper
	
	# define anat scan
	anat_scan=$anat_dir/"${subject}"_rec-NORM_T1w.nii.gz
	anat_scan=$anat_dir/anatUAC."${subject}".nii # this is the uniformized, intensity-corrected and ceiling-capped image.
	fs_anat=T1.nii

	# create anat folder in derivatives
	anat_deriv=$DIR/derivatives/"${subject}"/anat
	mkdir $anat_deriv

	# cd into subject folder
	cd $DIR/derivatives/FreeSurfer/"${subject}"/

	##### use SUMA to convert FreeSurfer output ##### 
	echo START WITH SUMA CONVERSION
	#@SUMA_Make_Spec_FS -sid ${subject} -NIFTI

	# cd into SUMA folder
	suma_dir=$DIR/derivatives/FreeSurfer/"${subject}"/SUMA
	cd $suma_dir

	##### create masks ##### 
	echo CREATE THE WM AND VENTRICLE MASK

	# define mask prefix
	GM_mask=FSmask_GM.nii.gz
	WM_mask=FSmask_WM.nii.gz
	vent_mask=FSmask_vent.nii.gz

	# define atlas
	desikan=aparc+aseg_REN_all.nii.gz
	destrieux=aparc.a2009s+aseg_REN_all.nii.gz
	destrieux_gm=aparc.a2009s+aseg_REN_gm.nii.gz

	# select ventricles from FS output aparc+aseg.nii (numbering corresponds to FS)
	3dcalc -a $destrieux -datum byte \
	-prefix $vent_mask -expr 'amongst(a,3,23)' 
	# "3" "Left-Lateral-Ventricle" 
	# "23" "Right-Lateral-Ventricle"

	# select the WM maps from FS output aparc+aseg.nii (numbering corresponds to FS)
	3dcalc -a $destrieux -datum byte \
	-prefix $WM_mask -expr 'amongst(a,1,5,21,25,41,42,43,44,45)'
	# "1" "Left-Cerebral-White-Matter" ""
	# "5" "Left-Cerebellum-White-Matter" ""
	# "21" "Right-Cerebral-White-Matter" ""
	# "25" "Right-Cerebellum-White-Matter" ""
	# "41" "CC_Posterior" ""
	# "42" "CC_Mid_Posterior" ""
	# "43" "CC_Central" ""
	# "44" "CC_Mid_Anterior" ""
	# "45" "CC_Anterior" ""
	# note: Chen et al. 2016 also included Brain-Stem, but this is omitted here

	# select GM mask
	3dmask_tool -input $destrieux_gm -prefix $GM_mask
	#3dmask_tool -input aparc+aseg_REN_wmat.nii.gz -prefix WM_mask.nii
	#3dmask_tool -input aparc+aseg_REN_vent.nii.gz -prefix VENT_mask.nii
done

