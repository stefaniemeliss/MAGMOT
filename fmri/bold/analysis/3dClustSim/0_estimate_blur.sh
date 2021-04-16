#!/bin/tcsh
source ~/.cshrc

################################################################################
# pre-processing for functional data using anfi_proc.py 
################################################################################

#module load afni19.3.03

# Set top level directory structure
set topdir = /storage/shared/research/cinn/2018/MAGMOT #study folder
echo $topdir

set derivroot = $topdir/derivatives
set outroot = $derivroot/afniproc

# define subject listecho $
set BIDSdir = $topdir/rawdata


cd $BIDSdir
set subjects	=(`ls -d sub*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-control*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-experimental04*`) # this creates an array containing all subjects in the BIDS directory
#echo $subjects


# determine variables to regulate flow
set run_task = true
set run_rest = false

#set subjects	= sub-control001

echo $#subjects

# set different smoothing kernels
set kernels = (8 0 4 6)
#set kernels = (8)

# determine clustsim directory
set clust_rest = $derivroot"/analysis/rest/3dClustSim"
mkdir $clust_rest
set clust_task = $derivroot"/analysis/magictrickwatching/3dClustSim"
mkdir $clust_rest

# for each subject in the subjects array
foreach subject ($subjects)

	echo $subject

	# Output directory: name for output
	set outdir  = $outroot/$subject
	cd $outdir # define PWD in which the script and results should get saved

    ###########################################################################
    ############################## RESTING-STATE ##############################
    ###########################################################################

    if ($run_rest == true) then

		set task = rest
		set subj = "$subject"_task-"$task"
	
		set output_dir = "$subj".results
		cd $output_dir
		set runs = (`count -digits 2 1 2`)

		foreach blur ( $kernels ) 

			set subjstr = "${subj}"_s"$blur"

			# ============================ blur estimation =============================
			# compute blur estimates
			touch blur_est.$subjstr.1D   # start with empty file

			# create directory for ACF curve files
			mkdir files_ACF_s"$blur"

			# -- estimate blur for each run in errts --
			touch blur.errts.$subjstr.1D

			# restrict to uncensored TRs, per run
			foreach run ( $runs )
				set trs = `1d_tool.py -infile X.s${blur}.xmat.1D -show_trs_uncensored encoded  \
						              -show_trs_run $run`
				if ( $trs == "" ) continue
				3dFWHMx -detrend -mask mask_epi_anat.$subj+tlrc                       \
						-ACF files_ACF_s"$blur"/out.3dFWHMx.ACF.errts.s$blur.r$run.1D                 \
						errts.$subjstr.fanaticor+tlrc"[$trs]" >> blur.errts.$subjstr.1D
			end

			# compute average FWHM blur (from every other row) and append
			set blurs = ( `3dTstat -mean -prefix - blur.errts.$subjstr.1D'{0..$(2)}'\'` )
			echo average errts FWHM blurs: $blurs
			echo "$blurs   # errts FWHM blur estimates" >> blur_est.$subjstr.1D

			# compute average ACF blur (from every other row) and append
			set blurs = ( `3dTstat -mean -prefix - blur.errts.$subjstr.1D'{1..$(2)}'\'` )
			echo average errts ACF blurs: $blurs
			echo "$blurs   # errts ACF blur estimates" >> blur_est.$subjstr.1D

			set out_file = out_FWHMx_s"$blur".1D

			# add data to output file
			if ($subject == sub-control001) then
				printf "$blurs" > "$clust_rest"/$out_file
			else
				printf "\n$blurs" >> "$clust_rest"/$out_file
			endif
		end

	endif

	cd $outdir

    ##################################################################################
    ############################## MAGIC TRICK WATCHING ##############################
    ##################################################################################

    if ($run_task == true) then

		set task = magictrickwatching
		set subj = "$subject"_task-"$task"
	
		set output_dir = "$subj".results
		cd $output_dir

		if ($subj == sub-experimental016) then
			set runs = (`count -digits 2 1 4`)
		else
			set runs = (`count -digits 2 1 3`)
		endif

		foreach blur ( $kernels ) 

			set subjstr = "${subj}"_s"$blur"

			# ============================ blur estimation =============================
			# compute blur estimates
			touch blur_est.$subjstr.1D   # start with empty file

			# create directory for ACF curve files
			mkdir files_ACF_s"$blur"

			# -- estimate blur for each run in errts --
			touch blur.errts.$subjstr.1D

			# restrict to uncensored TRs, per run
			foreach run ( $runs )
				set trs = `1d_tool.py -infile X.s${blur}.xmat.1D -show_trs_uncensored encoded  \
						              -show_trs_run $run`
				if ( $trs == "" ) continue
				3dFWHMx -detrend -mask mask_epi_anat.$subj+tlrc                       \
						-ACF files_ACF_s"$blur"/out.3dFWHMx.ACF.errts.s$blur.r$run.1D                 \
						errts.$subjstr.fanaticor+tlrc"[$trs]" >> blur.errts.$subjstr.1D
			end

			# compute average FWHM blur (from every other row) and append
			set blurs = ( `3dTstat -mean -prefix - blur.errts.$subjstr.1D'{0..$(2)}'\'` )
			echo average errts FWHM blurs: $blurs
			echo "$blurs   # errts FWHM blur estimates" >> blur_est.$subjstr.1D

			# compute average ACF blur (from every other row) and append
			set blurs = ( `3dTstat -mean -prefix - blur.errts.$subjstr.1D'{1..$(2)}'\'` )
			echo average errts ACF blurs: $blurs
			echo "$blurs   # errts ACF blur estimates" >> blur_est.$subjstr.1D

			set out_file = out_FWHMx_s"$blur".1D

			# add data to output file
			if ($subject == sub-control001) then
				printf "$blurs" > "$clust_task"/$out_file
			else
				printf "\n$blurs" >> "$clust_task"/$out_file
			endif

		end

	endif


end
