#!/bin/tcsh
source ~/.cshrc

################################################################################
# pre-processing for functional data using anfi_proc.py 
################################################################################

#module load afni19.3.03

# --------------------------------------------------------------------
# Script: s.2016_ChenEtal_02_ap.tcsh
#
# From:
# Chen GC, Taylor PA, Shin Y-W, Reynolds RC, Cox RW (2016). Untangling
# the Relatedness among Correlations, Part II: Inter-Subject
# Correlation Group Analysis through Linear Mixed-Effects
# Modeling. Neuroimage (in press).
#
# Originally run using: AFNI_16.1.16
# --------------------------------------------------------------------

# further modified based on Example 11. (see afni_proc.py -help)

# FMRI processing script, ISC movie data.
# Assumes previously run FS and SUMA commands, respectively: 
# $ recon-all -all -subject $subj -i $anat
# $ @SUMA_Make_Spec_FS -sid $subj -NIFTI

# Set top level directory structure
set topdir = /storage/shared/research/cinn/2018/MAGMOT #study folder
echo $topdir

set task = magictrickwatching
#set tasks = (magictrickwatching rest)

set derivroot = $topdir/derivatives
set fsroot = $derivroot/FreeSurfers
set outroot = $derivroot/afniproc

# define subject listecho $
set BIDSdir = $topdir/rawdata

cd $BIDSdir
set subjects	=(`ls -d sub*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-control*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-experimental*`) # this creates an array containing all subjects in the BIDS directory
#echo $subjects
echo $#subjects

#set subjects	= sub-experimental005
#set subjects	= sub-control002
#set subjects	= sub-control017



# for each subject in the subjects array
foreach subj ($subjects)

	#set subj	= "sub-experimental005"
	echo $subj

	# Output directory: name for output
	set outdir  = $outroot/$subj
	cd $outdir # define PWD in which the script and results should get saved

	# Input directory: unprocessed FMRI data
	set indir   = $BIDSdir/$subj/func

	# Input directory: anatomical derivatives
	set derivindir = $derivroot/$subj/anat
	set derivoutdir = $derivroot/$subj/func

	# Input directory: SSWarper output (anatomical non-linear aligned to MNI)
	set sswindir = $outdir/SSwarper

	# Input data: FreeSurfer results (anatomy, ventricle and WM maps)
	# all these files are in the ../derivatives/FreeSurfer/$SUBJ_ID/SUMA directory
	set anatSS = anatSS.${subj}.nii	
    set anatUAC = anatUAC.${subj}.nii

	set fsvent = ${subj}_space-orig_res-anat_label-vent_mask.nii.gz
	set fswm   = ${subj}_space-orig_res-anat_label-WM_mask.nii.gz
	set fsgm   = ${subj}_space-orig_res-anat_label-GM_mask.nii.gz
    set fsparc = ${subj}_space-orig_desc-Destrieux_dseg.nii.gz

    ####################################################
    ############### MAGIC TRICK WATCHING ###############
    ####################################################

    set task = magictrickwatching
    set subjstr = "$subj"_task-"$task"
    set POLORT = 6

	# Input data: list of partitioned EPIs (movie clips)
	set epi_dpattern = $outdir"/epi_b0_correct/"${subj}"_task-"${task}"_run-*_desc-b0corrected_EPI.nii.gz"
	echo $epi_dpattern

	if ($subj == sub-experimental016) then
		set epi_dpattern = ($indir"/sub-experimental016_task-magictrickwatching_run-1_desc-b0corrected_EPI.nii.gz" 		\
							$indir"/sub-experimental016_task-magictrickwatching_acq-1_run-2_desc-b0corrected_EPI.nii.gz" 	\
							$indir"/sub-experimental016_task-magictrickwatching_acq-2_run-2_desc-b0corrected_EPI.nii.gz"	\
							$indir"/sub-experimental016_task-magictrickwatching_run-3_desc-b0corrected_EPI.nii.gz")
	endif

    # determine minimal outlier in first run to use in volreg
    set first_run = $outdir"/epi_b0_correct/"${subj}"_task-"${task}"_run-1_desc-b0corrected_EPI.nii.gz"
    # run outlier count
    3dToutcount -automask -fraction -polort $POLORT -legendre  $first_run > outcount.first_run.1D
    # determine vol that has min outlier
    set min_out_first_run = `3dTstat -argmin -prefix - outcount.first_run.1D\'`
    # delete file
    rm outcount.first_run.1D

    # determine minimal outlier in first run to use in volreg
    set first_run = $indir"/"${subj}"_task-"$task"_run-1_bold.nii.gz"
    # run outlier count
    3dToutcount -automask -fraction -polort $POLORT -legendre  $first_run > outcount.first_run.1D
    # determine vol that has min outlier
    set min_out_first_run = `3dTstat -argmin -prefix - outcount.first_run.1D\'`
    # delete file
    rm outcount.first_run.1D

	# specify actual afni_proc.py
	afni_proc.py -subj_id "${subjstr}"					\
	    -blocks despike tshift align tlrc volreg mask blur scale regress          \
        -outlier_polort $POLORT                                             \
		-radial_correlate_blocks tcat volreg								\
	    -copy_anat $sswindir/$anatSS                                       \
		-anat_has_skull no													\
		-anat_follower anat_w_skull anat $sswindir/$anatUAC                 \
	    -anat_follower_ROI aaseg  anat $derivindir/$fsparc				\
	    -anat_follower_ROI aeseg  epi  $derivindir/$fsparc				\
	    -anat_follower_ROI FSvent epi  $derivindir/$fsvent                     \
	    -anat_follower_ROI FSWMe  epi  $derivindir/$fswm                       \
	    -anat_follower_ROI FSGMe  epi  $derivindir/$fsgm                       \
	    -anat_follower_erode FSvent FSWMe                                \
	    -dsets $epi_dpattern                                                \
	    -tcat_remove_first_trs 0                                            \
		-tshift_opts_ts -tpattern altplus									\
        -align_opts_aea -cost lpc+ZZ -giant_move -check_flip				\
	    -tlrc_base MNI152_2009_template_SSW.nii.gz		             		\
	    -tlrc_NL_warp                                                       \
		-tlrc_NL_warped_dsets $sswindir/anatQQ.$subj.nii 					\
			$sswindir/anatQQ.$subj.aff12.1D									\
			$sswindir/anatQQ."$subj"_WARP.nii								\
	    -volreg_base_ind 1 $min_out_first_run                               \
        -volreg_post_vr_allin yes                                           \
        -volreg_pvra_base_index MIN_OUTLIER                                 \
	    -volreg_align_e2a                                                   \
	    -volreg_tlrc_warp                                                   \
        -volreg_no_extent_mask                                             \
        -mask_dilate 8                                                      \
		-mask_epi_anat yes													\
		-blur_to_fwhm -blur_size 8                                          \
		-regress_motion_per_run												\
	    -regress_ROI_PC FSvent 3                                            \
		-regress_ROI_PC_per_run FSvent										\
	    -regress_make_corr_vols aeseg FSvent                                \
	    -regress_anaticor_fast                                              \
	    -regress_anaticor_label FSWMe                                       \
	    -regress_censor_motion 0.3                                          \
	    -regress_censor_outliers 0.1                                        \
	    -regress_apply_mot_types demean deriv                               \
	    -regress_est_blur_epits                                             \
	    -regress_est_blur_errts                                             \
	    -regress_run_clustsim no											\
		-regress_polort $POLORT												\
		-html_review_style pythonic

	# note: regress_polort 6 due to: 1140 volumns in total * 2 (TR = 2) = 2280 seconds
	# 2280s / 3 (number of runs) = 760s (run_length) per run [this matches with average block_durInSecs]
	# DEGREE: 1 + floor(run_length / 150.0)	= 1 + floor(5.066) = 6

	# trpattern = altplus due to slice_code = 3 (i.e. interleaved ascending) and slice_start = 0 in DICOM header

	# execute script
	tcsh -xef proc."${subjstr}" |& tee output.proc."${subjstr}"

    ##### compute subject mask that includes dilated GM (slightly touching WM & ventricle) in epi grid ###
    set output_dir = $subjstr.results

    # dilate subject anatomy mask (mask_anat.$subj_$task+tlrc is already in MNI space and on EPI grid)
    3dmask_tool -dilate_inputs 1 -inputs $output_dir/mask_anat."$subjstr"+tlrc. -prefix $output_dir/mask_anat_dilated."$subjstr"+tlrc.

    # remove eroded WM and eroded ventricle masks (anatomical follower, already in MI space & anatomical grid)
    3dcalc -a $output_dir/mask_anat_dilated."$subjstr"+tlrc. -b $output_dir/follow_ROI_FSWMe+tlrc. -c $output_dir/follow_ROI_FSvent+tlrc. -expr 'a-b-c' \
        -prefix $derivoutdir/"$subjstr"_space-MNI152NLin2009cAsym_res-epi_label-dilatedGM_mask.nii.gz

    # remove dilated subject anatomy mask
    rm $output_dir/mask_anat_dilated."$subjstr"+tlrc.*


    #############################################
    ############### RESTING-STATE ###############
    #############################################

    set task = rest
    set subjstr = "$subj"_task-"$task"
    set POLORT = 5

    # Input data: list of partitioned EPIs (resting state)
	set epi_dpattern = $outdir"/epi_b0_correct/"${subj}"_task-"${task}"_run-*_desc-b0corrected_EPI.nii.gz"
	echo $epi_dpattern

    # determine minimal outlier in first run to use in volreg
    set first_run = $outdir"/epi_b0_correct/"${subj}"_task-"${task}"_run-1_desc-b0corrected_EPI.nii.gz"
    # run outlier count
    3dToutcount -automask -fraction -polort $POLORT -legendre  $first_run > outcount.first_run.1D
    # determine vol that has min outlier
    set min_out_first_run = `3dTstat -argmin -prefix - outcount.first_run.1D\'`
    # delete file
    rm outcount.first_run.1D

	# specify actual afni_proc.py
	afni_proc.py -subj_id "${subjstr}"					\
	    -blocks despike tshift align tlrc volreg mask blur scale regress          \
        -outlier_polort $POLORT                                             \
		-radial_correlate_blocks tcat volreg								\
	    -copy_anat $sswindir/$anatSS                                       \
		-anat_has_skull no													\
		-anat_follower anat_w_skull anat $sswindir/$anatUAC                 \
	    -anat_follower_ROI aaseg  anat $derivindir/$fsparc				\
	    -anat_follower_ROI aeseg  epi  $derivindir/$fsparc				\
	    -anat_follower_ROI FSvent epi  $derivindir/$fsvent                     \
	    -anat_follower_ROI FSWMe  epi  $derivindir/$fswm                       \
	    -anat_follower_ROI FSGMe  epi  $derivindir/$fsgm                       \
	    -anat_follower_erode FSvent FSWMe                                \
	    -dsets $epi_dpattern                                                \
	    -tcat_remove_first_trs 0                                            \
		-tshift_opts_ts -tpattern altplus									\
        -align_opts_aea -cost lpc+ZZ -giant_move -check_flip				\
	    -tlrc_base MNI152_2009_template_SSW.nii.gz		             		\
	    -tlrc_NL_warp                                                       \
		-tlrc_NL_warped_dsets $sswindir/anatQQ.$subj.nii 					\
			$sswindir/anatQQ.$subj.aff12.1D									\
			$sswindir/anatQQ."$subj"_WARP.nii								\
	    -volreg_base_ind 1 $min_out_first_run                               \
        -volreg_post_vr_allin yes                                           \
        -volreg_pvra_base_index MIN_OUTLIER                                 \
	    -volreg_align_e2a                                                   \
	    -volreg_tlrc_warp                                                   \
        -volreg_no_extent_mask                                             \
        -mask_dilate 8                                                      \
		-mask_epi_anat yes													\
		-blur_to_fwhm -blur_size 8                                          \
		-regress_motion_per_run												\
	    -regress_ROI_PC FSvent 3                                            \
		-regress_ROI_PC_per_run FSvent										\
	    -regress_make_corr_vols aeseg FSvent                                \
	    -regress_anaticor_fast                                              \
	    -regress_anaticor_label FSWMe                                       \
	    -regress_censor_motion 0.3                                          \
	    -regress_censor_outliers 0.1                                        \
	    -regress_apply_mot_types demean deriv                               \
	    -regress_est_blur_epits                                             \
	    -regress_est_blur_errts                                             \
	    -regress_run_clustsim no											\
		-regress_polort $POLORT												\
		-html_review_style pythonic

	# note: regress_polort 5 due to: 600 volumns in total * 2 (TR = 2) = 1200 seconds
	# 1200 / 2 (number of runs) = 600s (run_length) per run
	# DEGREE: 1 + floor(run_length / 150.0)	= 1 + floor(4) = 5

	# trpattern = altplus due to slice_code = 3 (i.e. interleaved ascending) and slice_start = 0 in DICOM header

	# execute script
	tcsh -xef proc."${subjstr}" |& tee output.proc."${subjstr}"

    ##### compute subject mask that includes dilated GM (slightly touching WM & ventricle) in epi grid ###
    set output_dir = $subjstr.results

    # dilate subject anatomy mask (mask_anat.$subj_$task+tlrc is already in MNI space and on EPI grid)
    3dmask_tool -dilate_inputs 1 -inputs $output_dir/mask_anat."$subjstr"+tlrc. -prefix $output_dir/mask_anat_dilated."$subjstr"+tlrc.

    # remove eroded WM and eroded ventricle masks (anatomical follower, already in MI space & anatomical grid)
    3dcalc -a $output_dir/mask_anat_dilated."$subjstr"+tlrc. -b $output_dir/follow_ROI_FSWMe+tlrc. -c $output_dir/follow_ROI_FSvent+tlrc. -expr 'a-b-c' \
        -prefix $derivoutdir/"$subjstr"_space-MNI152NLin2009cAsym_res-epi_label-dilatedGM_mask.nii.gz

    # remove dilated subject anatomy mask
    rm $output_dir/mask_anat_dilated."$subjstr"+tlrc.*
							
end
