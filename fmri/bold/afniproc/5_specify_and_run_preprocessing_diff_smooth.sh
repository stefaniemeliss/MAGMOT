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

set derivroot = $topdir/derivatives
set fsroot = $derivroot/FreeSurfers
set outroot = $derivroot/afniproc

# define subject listecho $
set BIDSdir = $topdir/rawdata

cd $BIDSdir
set subjects	=(`ls -d sub*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-control*`) # this creates an array containing all subjects in the BIDS directory
set subjects	=(`ls -d sub-experimental04*`) # this creates an array containing all subjects in the BIDS directory
#echo $subjects
echo $#subjects

# determine variables to regulate flow
set run_preproc = true
set run_after = false

set subjects	= (sub-experimental024) 
set subjects	= sub-control002
#set subjects	= sub-control017

# set different smoothing kernels
set kernels = (4 6 8)
set kernels = (8)

# for each subject in the subjects array
foreach subject ($subjects)

	#set subject	= "sub-experimental005"
	echo $subject

	# Output directory: name for output
	set outdir  = $outroot/$subject
	cd $outdir # define PWD in which the script and results should get saved

	# Input directory: unprocessed FMRI data
	set indir   = $BIDSdir/$subject/func

	# Input directory: anatomical derivatives
	set derivindir = $derivroot/$subject/anat
	set derivoutdir = $derivroot/$subject/func
	set regressdir = $derivroot/$subject/regressors

	# Input directory: SSWarper output (anatomical non-linear aligned to MNI)
	set sswindir = $outdir/SSwarper

	# Input data: FreeSurfer results (anatomy, ventricle and WM maps)
	# all these files are in the ../derivatives/FreeSurfer/$SUBJ_ID/SUMA directory
	#set anatSS = anatSS.${subject}.nii	
	set anatQQ = anatQQ.${subject}.nii	
    #set anatUAC = anatUAC.${subject}.nii

	#set fsvent = ${subject}_space-orig_res-anat_label-vent_mask.nii.gz
	#set fswm   = ${subject}_space-orig_res-anat_label-WM_mask.nii.gz
	#set fsgm   = ${subject}_space-orig_res-anat_label-GM_mask.nii.gz
    #set fsparc = ${subject}_space-orig_desc-Destrieux_dseg.nii.gz

    ##################################################################################
    ############################## MAGIC TRICK WATCHING ##############################
    ##################################################################################

    set task = magictrickwatching
    set subjstr = "$subject"_task-"$task"
    set POLORT = 6

	############################# SET UP PRE-PROCESSING ############################# 

    if ($run_preproc == true) then

	    # Input data: list of partitioned EPIs (movie clips)
	    set epi_dpattern = ${subjstr}".results/pb03."${subjstr}".r0*.volreg+tlrc.HEAD"

	    echo $epi_dpattern

		set blur = 0

		set subj = "${subjstr}"_s"$blur"

	    # specify actual afni_proc.py
	    afni_proc.py -subj_id $subj											\
	        -blocks mask scale regress          								\
			-copy_files $regressdir/"$subjstr"_run-1_label-ventriclePC_regressor.1D	\
					$regressdir/"$subjstr"_run-2_label-ventriclePC_regressor.1D	\
					$regressdir/"$subjstr"_run-3_label-ventriclePC_regressor.1D	\
					$regressdir/"$subjstr"_label-localWM_regressor.nii.gz		\
			-outlier_count no													\
	        -copy_anat $sswindir/$anatQQ                                       	\
		    -anat_has_skull no													\
	        -dsets $epi_dpattern                                                \
            -mask_dilate 8                                                      \
		    -mask_epi_anat yes													\
		    -regress_motion_file $regressdir/"$subjstr"_label-mot_regressor.1D	\
		    -regress_motion_per_run												\
	        -regress_censor_extern $regressdir/"$subjstr"_label-censorTRs_regressor.1D		\
	        -regress_apply_mot_types demean deriv                               \
			-regress_opts_3dD -ortvec "$subjstr"_run-1_label-ventriclePC_regressor.1D	ventPC_r01	\
					-ortvec "$subjstr"_run-2_label-ventriclePC_regressor.1D	ventPC_r02	\
					-ortvec "$subjstr"_run-3_label-ventriclePC_regressor.1D	ventPC_r03	\
	        -regress_run_clustsim no											\
		    -regress_polort $POLORT												\
			-regress_compute_tsnr no											\
			-regress_compute_gcor no											\
			-no_epi_review

		# regress_censor_extern --> combined motion & outlier censor

	    # execute script
	    #tcsh -xef proc."${subj}" |& tee output.proc."${subj}"

		#set output_dir = "$subj".results
		#cd $output_dir
		#set runs = (`count -digits 2 1 3`)

		# note TRs that were not censored
		#set ktrs = `1d_tool.py -infile "$subjstr"_label-censorTRs_regressor.1D		\
        #               -show_trs_uncensored encoded`


		# ============================ local WM regres =============================
		# -- use 3dTproject to project out regression matrix --
		#    (make errts like 3dDeconvolve, but more quickly)
		#3dTproject -polort 0 -input pb0*.$subj.r*.scale+tlrc.HEAD					\
		#		   -censor "$subjstr"_label-censorTRs_regressor.1D -cenmode ZERO	\
		#		   -dsort "$subjstr"_label-localWM_regressor.nii.gz					\
		#		   -ort X.nocensor.xmat.1D -prefix errts.$subj.fanaticor

		# ============================ blur estimation =============================
		# compute blur estimates
		#touch blur_est.$subj.1D   # start with empty file

		# create directory for ACF curve files
		#mkdir files_ACF

		# -- estimate blur for each run in epits --
		#touch blur.epits.1D

		# restrict to uncensored TRs, per run
		#foreach run ( $runs )
		#	set trs = `1d_tool.py -infile X.xmat.1D -show_trs_uncensored encoded  \
		#		                  -show_trs_run $run`
		#	if ( $trs == "" ) continue
		#	3dFWHMx -detrend -mask mask_epi_anat.$subj+tlrc                       \
		#		    -ACF files_ACF/out.3dFWHMx.ACF.epits.r$run.1D                 \
		#		    all_runs.$subj+tlrc"[$trs]" >> blur.epits.1D
		#end

		# compute average FWHM blur (from every other row) and append
		#set blurs = ( `3dTstat -mean -prefix - blur.epits.1D'{0..$(2)}'\'` )
		#echo average epits FWHM blurs: $blurs
		#echo "$blurs   # epits FWHM blur estimates" >> blur_est.$subj.1D

		# compute average ACF blur (from every other row) and append
		#set blurs = ( `3dTstat -mean -prefix - blur.epits.1D'{1..$(2)}'\'` )
		#echo average epits ACF blurs: $blurs
		#echo "$blurs   # epits ACF blur estimates" >> blur_est.$subj.1D

		# -- estimate blur for each run in errts --
		#touch blur.errts.1D

		# restrict to uncensored TRs, per run
		#foreach run ( $runs )
		#	set trs = `1d_tool.py -infile X.xmat.1D -show_trs_uncensored encoded  \
		#		                  -show_trs_run $run`
		#	if ( $trs == "" ) continue
		#	3dFWHMx -detrend -mask mask_epi_anat.$subj+tlrc                       \
		#		    -ACF files_ACF/out.3dFWHMx.ACF.errts.r$run.1D                 \
		#		    errts.$subj.fanaticor+tlrc"[$trs]" >> blur.errts.1D
		#end

		# compute average FWHM blur (from every other row) and append
		#set blurs = ( `3dTstat -mean -prefix - blur.errts.1D'{0..$(2)}'\'` )
		#echo average errts FWHM blurs: $blurs
		#echo "$blurs   # errts FWHM blur estimates" >> blur_est.$subj.1D

		# compute average ACF blur (from every other row) and append
		#set blurs = ( `3dTstat -mean -prefix - blur.errts.1D'{1..$(2)}'\'` )
		#echo average errts ACF blurs: $blurs
		#echo "$blurs   # errts ACF blur estimates" >> blur_est.$subj.1D

	    # define files: final preprocessing output is errts.$subjects.magictrickwatching_perRun.fanaticor
        #set file_string = errts."$subj".fanaticor+tlrc.HEAD

	    #set file = (`ls $file_string | xargs -n 1 basename`)

	    # define niifile
	    #set niifile = "$subjstr""_desc-s""$blur"preproc_bold.nii.gz

	    # do the AFNI to .nii conversion
	    #3dAFNItoNIFTI -prefix $derivoutdir/$niifile $output_dir/$file_string

		#cd $outdir

    foreach blur ( $kernels )

		set subj = "${subjstr}"_s"$blur"

	    # specify actual afni_proc.py
	    afni_proc.py -subj_id $subj											\
	        -blocks mask blur scale regress          								\
			-copy_files $regressdir/"$subjstr"_run-1_label-ventriclePC_regressor.1D	\
					$regressdir/"$subjstr"_run-2_label-ventriclePC_regressor.1D	\
					$regressdir/"$subjstr"_run-3_label-ventriclePC_regressor.1D	\
					$regressdir/"$subjstr"_label-localWM_regressor.nii.gz		\
			-outlier_count no													\
	        -copy_anat $sswindir/$anatQQ                                       	\
		    -anat_has_skull no													\
	        -dsets $epi_dpattern                                                \
            -mask_dilate 8                                                      \
		    -mask_epi_anat yes													\
            -blur_to_fwhm -blur_size $blur											\
		    -regress_motion_file $regressdir/"$subjstr"_label-mot_regressor.1D	\
		    -regress_motion_per_run												\
	        -regress_censor_extern $regressdir/"$subjstr"_label-censorTRs_regressor.1D		\
			-regress_opts_3dD -ortvec "$subjstr"_run-1_label-ventriclePC_regressor.1D	ventPC_r01	\
					-ortvec "$subjstr"_run-2_label-ventriclePC_regressor.1D	ventPC_r02	\
					-ortvec "$subjstr"_run-3_label-ventriclePC_regressor.1D	ventPC_r03	\
	        -regress_apply_mot_types demean deriv                               \
	        -regress_run_clustsim no											\
		    -regress_polort $POLORT												\
			-regress_compute_tsnr no											\
			-regress_compute_gcor no											\
			-no_epi_review


		# regress_censor_extern --> combined motion & outlier censor

	    # execute script
	    tcsh -xef proc."${subj}" |& tee output.proc."${subj}"

		set output_dir = "$subj".results
		cd $output_dir
		set runs = (`count -digits 2 1 3`)

		# note TRs that were not censored
		set ktrs = `1d_tool.py -infile "$subjstr"_label-censorTRs_regressor.1D		\
                       -show_trs_uncensored encoded`


		# ============================ local WM regres =============================
		# -- use 3dTproject to project out regression matrix --
		#    (make errts like 3dDeconvolve, but more quickly)
		3dTproject -polort 0 -input pb0*.$subj.r*.scale+tlrc.HEAD					\
				   -censor "$subjstr"_label-censorTRs_regressor.1D -cenmode ZERO	\
				   -dsort "$subjstr"_label-localWM_regressor.nii.gz					\
				   -ort X.nocensor.xmat.1D -prefix errts.$subj.fanaticor

		# ============================ blur estimation =============================
		# compute blur estimates
		touch blur_est.$subj.1D   # start with empty file

		# create directory for ACF curve files
		mkdir files_ACF

		# -- estimate blur for each run in epits --
		touch blur.epits.1D

		# restrict to uncensored TRs, per run
		foreach run ( $runs )
			set trs = `1d_tool.py -infile X.xmat.1D -show_trs_uncensored encoded  \
				                  -show_trs_run $run`
			if ( $trs == "" ) continue
			3dFWHMx -detrend -mask mask_epi_anat.$subj+tlrc                       \
				    -ACF files_ACF/out.3dFWHMx.ACF.epits.r$run.1D                 \
				    all_runs.$subj+tlrc"[$trs]" >> blur.epits.1D
		end

		# compute average FWHM blur (from every other row) and append
		set blurs = ( `3dTstat -mean -prefix - blur.epits.1D'{0..$(2)}'\'` )
		echo average epits FWHM blurs: $blurs
		echo "$blurs   # epits FWHM blur estimates" >> blur_est.$subj.1D

		# compute average ACF blur (from every other row) and append
		set blurs = ( `3dTstat -mean -prefix - blur.epits.1D'{1..$(2)}'\'` )
		echo average epits ACF blurs: $blurs
		echo "$blurs   # epits ACF blur estimates" >> blur_est.$subj.1D

		# -- estimate blur for each run in errts --
		touch blur.errts.1D

		# restrict to uncensored TRs, per run
		foreach run ( $runs )
			set trs = `1d_tool.py -infile X.xmat.1D -show_trs_uncensored encoded  \
				                  -show_trs_run $run`
			if ( $trs == "" ) continue
			3dFWHMx -detrend -mask mask_epi_anat.$subj+tlrc                       \
				    -ACF files_ACF/out.3dFWHMx.ACF.errts.r$run.1D                 \
				    errts.$subj.fanaticor+tlrc"[$trs]" >> blur.errts.1D
		end

		# compute average FWHM blur (from every other row) and append
		set blurs = ( `3dTstat -mean -prefix - blur.errts.1D'{0..$(2)}'\'` )
		echo average errts FWHM blurs: $blurs
		echo "$blurs   # errts FWHM blur estimates" >> blur_est.$subj.1D

		# compute average ACF blur (from every other row) and append
		set blurs = ( `3dTstat -mean -prefix - blur.errts.1D'{1..$(2)}'\'` )
		echo average errts ACF blurs: $blurs
		echo "$blurs   # errts ACF blur estimates" >> blur_est.$subj.1D

	    ##### convert errts.anaticor output to nifti and copy to derivatives directory #####

		# define files: final preprocessing output is errts.$subjects.magictrickwatching_perRun.fanaticor
        set file_string = errts."$subj".fanaticor+tlrc.HEAD

	    #set file = (`ls $file_string | xargs -n 1 basename`)

	    # define niifile
	    set niifile = "$subjstr""_desc-s""$blur"preproc_bold.nii.gz

	    # do the AFNI to .nii conversion
	    3dAFNItoNIFTI -prefix $derivoutdir/$niifile $output_dir/$file_string

		cd $outdir

    endif

#** ERROR: (FAILED) attempt to over-write file 3dFWHMx.1D

    ###########################################################################
    ############################## RESTING-STATE ##############################
    ###########################################################################

    #set task = rest
    #set subjstr = "$subj"_task-"$task"
    #set POLORT = 5

							
end
