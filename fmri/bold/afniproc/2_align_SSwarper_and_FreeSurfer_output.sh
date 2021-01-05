#!/bin/bash

################################################################################
# converting FreeSurfer to AFNI to create masks for pre-processing
################################################################################

# define DIR
DIR=/storage/shared/research/cinn/2018/MAGMOT

# change directory to the raw NIFI files
cd $DIR/MAGMOT_BIDS/

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
#subjects=(sub-experimental005)
#subjects=(sub-control001)

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

	# define SSwarper output folder
	anat_dir=$DIR/derivatives/afniproc/"${subject}"/SSwarper
	
	# define anat scan
	anat_scan=$anat_dir/anatUAC."${subject}".nii # this is the uniformized, intensity-corrected and ceiling-capped image.
/storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/sub-control001/SSwarper/anatUAC.sub-control001.nii
	fs_anat=T1.nii

	# create anat folder in derivatives
	anat_deriv=$DIR/derivatives/"${subject}"/anat
	mkdir $anat_deriv

	# cd into SUMA folder
	suma_dir=$DIR/derivatives/FreeSurfer/"${subject}"/SUMA
	cd $suma_dir

	##### align SUMA and anatomy ######
	echo ALIGN T1w and SUMA CONVERSION OUTPUT

	# define mask prefix
	GM_mask=FSmask_GM.nii.gz
	WM_mask=FSmask_WM.nii.gz
	vent_mask=FSmask_vent.nii.gz
	# define atlas
	desikan=aparc+aseg_REN_all.nii.gz
	destrieux=aparc.a2009s+aseg_REN_all.nii.gz

	# problem: the original T1w scan () and the output of the conversion (e.g. T1 etc) are not aligned
    # for the WM and ventricle masks to be used within afni_proc.py (using -anat_follower), they MUST be aligned with the original anatomy input (i.e. anatSS.$subject.nii)

    # 1. align center of FS T1 to UAC.anat (and apply same parameters to the masks) --> this creates *_shft+orig.
    @Align_Centers  \
        -cm \
        -dset $fs_anat \
        -base $anat_scan \
		-child $GM_mask $WM_mask $vent_mask $destrieux $desikan

	# 2. align shifted FS T1 to anatomical scan (and apply same parameters to the masks) --> this creates *_shft+orig.
	fs_shft=T1_shft.nii
	wm_shft=FSmask_WM_shft.nii.gz
	gm_shft=FSmask_GM_shft.nii.gz
	vent_shft=FSmask_vent_shft.nii.gz
	desikan_shft=aparc+aseg_REN_all_shft.nii.gz
	destrieux_shft=aparc.a2009s+aseg_REN_all_shft.nii.gz

	align_epi_anat.py			\
		-dset1 $anat_scan		\
		-dset2 $fs_shft         \
		-dset2to1				\
		-giant_move				\
		-Allineate_opts '-final NN' \
		-child_dset2 $gm_shft $wm_shft $vent_shft $desikan_shft $destrieux_shft

	# 4. check alignment by overlaying anatomical scan used in SSwarper (anatUAC) and masks
	# copy anatomical scan
	anatUAC=anatUAC.nii.gz
	3dcopy $anat_scan ./$anatUAC 

	# define file names
	wm_al=FSmask_WM_shft_al+orig.
	gm_al=FSmask_GM_shft_al+orig.
	vent_al=FSmask_vent_shft_al+orig.

	#fs_al=T1_shft_al+orig.
	desikan_al=aparc+aseg_REN_all_shft_al+orig.
	destrieux_al=aparc.a2009s+aseg_REN_all_shft_al+orig.

	# define output file names
	out_gm=overlay.$subject.gm
	out_wm=overlay.$subject.wm
	out_vent=overlay.$subject.vent

	# create image overlaying T1w and grey matter mask
    afni -noplugins -no_detach                               \
         -com "OPEN_WINDOW sagittalimage opacity=4"          \
         -com "OPEN_WINDOW axialimage opacity=4"             \
         -com "OPEN_WINDOW coronalimage opacity=4"           \
         -com "SWITCH_UNDERLAY $anatUAC"                          \
         -com "SWITCH_OVERLAY $gm_al"                           \
         -com "SEE_OVERLAY +"                                \
         -com "SAVE_JPEG sagittalimage $out_gm.sag.jpg blowup=4"     \
         -com "SAVE_JPEG coronalimage  $out_gm.cor.jpg blowup=4"     \
         -com "SAVE_JPEG axialimage    $out_gm.axi.jpg blowup=4"     \
         -com "QUITT"                                        \
       $suma_dir


	# create image overlaying T1w and white matter mask
    afni -noplugins -no_detach                               \
         -com "OPEN_WINDOW sagittalimage opacity=4"          \
         -com "OPEN_WINDOW axialimage opacity=4"             \
         -com "OPEN_WINDOW coronalimage opacity=4"           \
         -com "SWITCH_UNDERLAY $anatUAC"                          \
         -com "SWITCH_OVERLAY $wm_al"                           \
         -com "SEE_OVERLAY +"                                \
         -com "SAVE_JPEG sagittalimage $out_wm.sag.jpg blowup=4"     \
         -com "SAVE_JPEG coronalimage  $out_wm.cor.jpg blowup=4"     \
         -com "SAVE_JPEG axialimage    $out_wm.axi.jpg blowup=4"     \
         -com "QUITT"                                        \
       $suma_dir

	# create image overlaying T1w and ventricle mask
    afni -noplugins -no_detach                               \
         -com "OPEN_WINDOW sagittalimage opacity=4"          \
         -com "OPEN_WINDOW axialimage opacity=4"             \
         -com "OPEN_WINDOW coronalimage opacity=4"           \
         -com "SWITCH_UNDERLAY $anatUAC"                          \
         -com "SWITCH_OVERLAY $vent_al"                           \
         -com "SEE_OVERLAY +"                                \
         -com "SAVE_JPEG sagittalimage $out_vent.sag.jpg blowup=4"     \
         -com "SAVE_JPEG coronalimage  $out_vent.cor.jpg blowup=4"     \
         -com "QUITT"                                        \
       $suma_dir


    # 4. copy aligned anatomical masks (and convert them .to nii.gz)
    3dcopy $gm_al $anat_deriv/"${subject}"_space-orig_res-anat_label-GM_mask.nii.gz
    3dcopy $wm_al $anat_deriv/"${subject}"_space-orig_res-anat_label-WM_mask.nii.gz
    3dcopy $vent_al $anat_deriv/"${subject}"_space-orig_res-anat_label-vent_mask.nii.gz

	# copy discrete segmentation files to derivatives/${subject}/anat folder
	# note: REN = renumbered
	# for Desikan-Killiany, see https://afni.nimh.nih.gov/pub/dist/src/scripts_install/afni_fs_aparc+aseg_2000.txt 
	# for Destrieux, see https://afni.nimh.nih.gov/pub/dist/src/scripts_install/afni_fs_aparc+aseg_2009.txt
	3dcopy $desikan_al $anat_deriv/"${subject}"_space-orig_desc-Desikan-Killiany_dseg.nii.gz
	3dcopy $destrieux_al $anat_deriv/"${subject}"_space-orig_desc-Destrieux_dseg.nii.gz

done










