#!/bin/tcsh
source ~/.cshrc

################################################################################
# pre-processing for functional data using anfi_proc.py 
################################################################################

module load afni19.3.03

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
set derivroot = $topdir/derivatives
set fsroot = $derivroot/FreeSurfers
set outroot = $derivroot/afniproc

# define subject listecho $
set BIDSdir = $topdir/MAGMOT_BIDS

cd $BIDSdir
set subjects	=(`ls -d sub*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-control*`) # this creates an array containing all subjects in the BIDS directory
set subjects	=(`ls -d sub-experimental*`) # this creates an array containing all subjects in the BIDS directory
echo $subjects
echo $#subjects

set subjects	= sub-experimental016
#set subjects	= sub-control001



# for each subject in the subjects array
foreach subj ($subjects)

	#set subj	= "sub-experimental005"
	echo $subj

	# Output directory: name for output
	set outdir  = $outroot/$subj
	cd $outdir # define PWD in which the script and results should get saved

	# Input directory: unprocessed FMRI data
	set indir   = $derivroot/$subj/func

	# Input directory: anatomical derivatives
	set derivindir = $derivroot/$subj/anat

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

	# Input data: list of partitioned EPIs (movie clips)
	set epi_dpattern = $indir"/"${subj}"_task-magictrickwatching_run-*_desc-cut_bold.nii.gz"
	echo $epi_dpattern

	if ($subj == sub-experimental016) then
		set epi_dpattern = ($indir"/sub-experimental016_task-magictrickwatching_run-1_desc-cut_bold.nii.gz" 		\
							$indir"/sub-experimental016_task-magictrickwatching_acq-1_run-2_desc-cut_bold.nii.gz" 	\
							$indir"/sub-experimental016_task-magictrickwatching_acq-2_run-2_desc-cut_bold.nii.gz"	\
							$indir"/sub-experimental016_task-magictrickwatching_run-3_desc-cut_bold.nii.gz")
	endif

echo 0000

	# specify actual afni_proc.py
	afni_proc.py -subj_id "${subj}"_task-"${task}"					\
	    -blocks despike tshift align tlrc volreg blur mask regress          \
		-radial_correlate_blocks tcat volreg								\
	    -copy_anat $sswindir/$anatSS                                       \
		-anat_has_skull no													\
		-anat_follower anat_w_skull anat $sswindir/$anatUAC                 \
	    -anat_follower_ROI aaseg  anat $derivindir/$fsparc				\
	    -anat_follower_ROI aeseg  epi  $derivindir/$fsparc				\
	    -anat_follower_ROI FSvent epi  $derivindir/$fsvent                     \
	    -anat_follower_ROI FSWMe  epi  $derivindir/$fswm                       \
	    -anat_follower_ROI FSGMe  epi  $derivindir/$fsgm                       \
	    -anat_follower_erode FSvent FSWMe FSGMe                             \
	    -dsets $epi_dpattern                                                \
	    -tcat_remove_first_trs 0                                            \
		-tshift_opts_ts -tpattern altplus									\
        -align_opts_aea -cost lpc+ZZ -giant_move -check_flip				\
	    -tlrc_base MNI152_2009_template_SSW.nii.gz		             		\
	    -tlrc_NL_warp                                                       \
		-tlrc_NL_warped_dsets $sswindir/anatQQ.$subj.nii 					\
			$sswindir/anatQQ.$subj.aff12.1D									\
			$sswindir/anatQQ."$subj"_WARP.nii								\
	    -volreg_align_to MIN_OUTLIER                                        \
        -volreg_post_vr_allin yes                                           \
        -volreg_pvra_base_index MIN_OUTLIER                                 \
	    -volreg_align_e2a                                                   \
	    -volreg_tlrc_warp                                                   \
		-blur_to_fwhm -blur_size 8 -blur_in_automask						\
		-mask_epi_anat yes													\
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
		-regress_polort 6													\
		-html_review_style pythonic

	# note: regress_polort 6 due to: 1140 volumns in total * 2 (TR = 2) = 2280 seconds
	# 2280s / 3 (number of runs) = 760s (run_length) per run [this matches with average block_durInSecs]
	# DEGREE: 1 + floor(run_length / 150.0)	= 1 + floor(5.066) = 6

	# trpattern = altplus due to slice_code = 3 (i.e. interleaved ascending) and slice_start = 0 in DICOM header

	# execute script
	tcsh -xef proc."${subj}"_task-"${task}" |& tee output.proc."${subj}"_task-"${task}"

							
end
#done
