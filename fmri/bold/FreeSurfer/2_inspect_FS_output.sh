#!/bin/bash

# load FreeSurfer
module load freesurfer6.0.0
source /usr/share/freesurfer/SetUpFreeSurfer.sh 

# change dir
export FREESURFER_HOME=/usr/share/freesurfer/
export SUBJECTS_DIR=/storage/shared/research/cinn/2018/MAGMOT/MAGMOT_BIDS

# change directory to the raw NIFTI files
cd /storage/shared/research/cinn/2018/MAGMOT/MAGMOT_BIDS/


# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
subjects=(sub-control035)

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	# define subject dir
	sub_dir=/storage/shared/research/cinn/2018/MAGMOT/derivatives/FreeSurfer/"${subject}"_edit/
	#sub_dir=/storage/shared/research/cinn/2018/MAGMOT/derivatives/FreeSurfer/sub-control003/
	cd $sub_dir

	# viewing volumes (-v) and surfaces (-f) with Freeview
	#freeview -v \
	#mri/brainmask.mgz \
	#mri/T1.mgz \
	#mri/wm.mgz:opacity=0.7 \
	#mri/aseg.mgz:colormap=lut:opacity=0.4 \
	#-f \
	#surf/lh.white:edgecolor=yellow \
	#surf/lh.pial:edgecolor=red \
	#surf/rh.white:edgecolor=yellow \
	#surf/rh.pial:edgecolor=red

	#surf/lh.white:edgecolor=blue \
	#surf/rh.white:edgecolor=blue \
	#surf/lh.pial:annot=aparc.a2009s.annot:name=pial_aparc_des:visible=1 \
	#surf/rh.pial:annot=aparc.a2009s.annot:name=pial_aparc_des:visible=1 \
	#--viewport 3d --hide-3d-slices

	freeview -v mri/brainmask.mgz \
	mri/wm.mgz:colormap=heat:opacity=0.4 \
	-f surf/lh.white:edgecolor=yellow \
	surf/lh.pial:edgecolor=red \
	surf/rh.white:edgecolor=yellow \
	surf/rh.pial:edgecolor=red \
	surf/rh.orig.nofix \
	surf/rh.orig:edgecolor=blue \
	surf/lh.orig.nofix \
	surf/lh.orig:edgecolor=blue 

done

