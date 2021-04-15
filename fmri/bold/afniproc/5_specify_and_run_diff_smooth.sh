#!/bin/tcsh
source ~/.cshrc

################################################################################
# pre-processing for functional data using anfi_proc.py 
################################################################################

module load afni19.3.03

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

set subjects	= (sub-experimental016) 
#set subjects	= sub-control002
#set subjects	= sub-control017

#set subjects	= (sub-experimental004 sub-experimental005 sub-experimental006 sub-experimental008 sub-experimental010 sub-experimental012 sub-experimental014 sub-experimental016 sub-experimental018 sub-experimental020 sub-experimental022 sub-experimental024 sub-experimental026 sub-experimental028 sub-experimental030)

#set subjects	= (sub-experimental032 sub-experimental034 sub-experimental036 sub-experimental038 sub-experimental040 sub-experimental042 sub-experimental044 sub-experimental046 sub-experimental048 sub-experimental050)
echo $#subjects

# set different smoothing kernels
set kernels = (4 6 8)
#set kernels = (8)

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

    ##################################################################################
    ############################## MAGIC TRICK WATCHING ##############################
    ##################################################################################

    set task = magictrickwatching
    set subj = "$subject"_task-"$task"
    set POLORT = 6
	
	set output_dir = "$subj".results
	cd $output_dir
	if ($subject == sub-experimental016) then
		set runs = (`count -digits 2 1 4`)
	else
		set runs = (`count -digits 2 1 3`)
	endif

	################################### NO BLUR	################################### 

	set blur = 0
	set subjstr = "${subj}"_s"$blur"

    if ($run_task == true) then

		# ================================= scale ==================================
		# scale each voxel time series to have a mean of 100
		# (be sure no negatives creep in)
		# (subject to a range of [0,200])
		foreach run ( $runs )
			3dTstat -prefix rm.mean_r$run pb03.$subj.r$run.volreg+tlrc
			3dcalc -a pb03.$subj.r$run.volreg+tlrc -b rm.mean_r$run+tlrc \
			       -c mask_epi_anat.$subj+tlrc                           \
			       -expr 'c * min(200, a/b*100)*step(a)*step(b)'         \
			       -prefix pb00.$subjstr.r$run.scale
		end

		# ================================ regress =================================
		# ------------------------------
		# run the regression analysis

	    if ($subject == sub-experimental016) then

			3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                         \
				-mask mask_epi_anat.$subj+tlrc                                        \
				-censor censor_${subj}_combined_2.1D                                  \
				-ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
				-ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
				-ortvec ROIPC.FSvent.r03.1D ROIPC.FSvent.r03                          \
				-ortvec ROIPC.FSvent.r04.1D ROIPC.FSvent.r04                          \
				-ortvec mot_demean.r01.1D mot_demean_r01                              \
				-ortvec mot_demean.r02.1D mot_demean_r02                              \
				-ortvec mot_demean.r03.1D mot_demean_r03                              \
				-ortvec mot_demean.r04.1D mot_demean_r04                              \
				-ortvec mot_deriv.r01.1D mot_deriv_r01                                \
				-ortvec mot_deriv.r02.1D mot_deriv_r02                                \
				-ortvec mot_deriv.r03.1D mot_deriv_r03                                \
				-ortvec mot_deriv.r04.1D mot_deriv_r04                                \
				-polort $POLORT                                                             \
				-num_stimts 0                                                         \
				-fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
				-x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
				-fitts fitts.$subjstr                                                 \
				-errts errts.${subjstr}                                               \
				-x1D_stop                                                             \
				-bucket stats.$subjstr

		else

			3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                         \
				-mask mask_epi_anat.$subj+tlrc                                        \
				-censor censor_${subj}_combined_2.1D                                  \
				-ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
				-ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
				-ortvec ROIPC.FSvent.r03.1D ROIPC.FSvent.r03                          \
				-ortvec mot_demean.r01.1D mot_demean_r01                              \
				-ortvec mot_demean.r02.1D mot_demean_r02                              \
				-ortvec mot_demean.r03.1D mot_demean_r03                              \
				-ortvec mot_deriv.r01.1D mot_deriv_r01                                \
				-ortvec mot_deriv.r02.1D mot_deriv_r02                                \
				-ortvec mot_deriv.r03.1D mot_deriv_r03                                \
				-polort $POLORT                                                             \
				-num_stimts 0                                                         \
				-fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
				-x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
				-fitts fitts.$subjstr                                                 \
				-errts errts.${subjstr}                                               \
				-x1D_stop                                                             \
				-bucket stats.$subjstr

		endif


		# -- use 3dTproject to project out regression matrix --
		#    (make errts like 3dDeconvolve, but more quickly)
		3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
				   -mask mask_epi_anat.$subj+tlrc                                 \
				   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
				   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.tproject

		# if 3dDeconvolve fails, terminate the script
		if ( $status != 0 ) then
			echo '---------------------------------------'
			echo '** 3dDeconvolve error, failing...'
			echo '   (consider the file 3dDeconvolve.err)'
			exit
		endif

		# -- use 3dTproject to project out regression matrix --
		#    (make errts like 3dDeconvolve, but more quickly)
		3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
				   -mask mask_epi_anat.$subj+tlrc                                 \
				   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
				   -dsort Local_FSWMe_rall+tlrc                                   \
				   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.fanaticor

		# ================================ convert =================================

		# define niifile
		set niifile = $subj"_desc-s"$blur"preproc_bold.nii.gz"

		# do the AFNI to .nii conversion
		3dAFNItoNIFTI -prefix $derivoutdir/$niifile errts.${subjstr}.fanaticor+tlrc

		# remove temporary files
		\rm -f rm.*

		###################################   BLUR	################################### 

		foreach blur ( $kernels ) 

			set subjstr = "${subj}"_s"$blur"

			# ================================== blur ==================================
			# blur each volume of each run
			foreach run ( $runs )
				3dBlurToFWHM -FWHM $blur -mask mask_epi_anat.$subj+tlrc \
						     -input pb03.$subj.r$run.volreg+tlrc    \
						     -prefix pb00.$subjstr.r$run.blur 
			end

			# ================================= scale ==================================
			# scale each voxel time series to have a mean of 100
			# (be sure no negatives creep in)
			# (subject to a range of [0,200])
			foreach run ( $runs )
				3dTstat -prefix rm.mean_r$run pb00.$subjstr.r$run.blur+tlrc
				3dcalc -a pb00.$subjstr.r$run.blur+tlrc -b rm.mean_r$run+tlrc \
					   -expr 'min(200, a/b*100)*step(a)*step(b)'           \
					   -prefix pb00.$subjstr.r$run.scale
			end

			# ================================ regress =================================
			# ------------------------------
			# run the regression analysis
			3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                         \
				-mask mask_epi_anat.$subj+tlrc                                        \
				-censor censor_${subj}_combined_2.1D                                  \
				-ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
				-ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
				-ortvec ROIPC.FSvent.r03.1D ROIPC.FSvent.r03                          \
				-ortvec mot_demean.r01.1D mot_demean_r01                              \
				-ortvec mot_demean.r02.1D mot_demean_r02                              \
				-ortvec mot_demean.r03.1D mot_demean_r03                              \
				-ortvec mot_deriv.r01.1D mot_deriv_r01                                \
				-ortvec mot_deriv.r02.1D mot_deriv_r02                                \
				-ortvec mot_deriv.r03.1D mot_deriv_r03                                \
				-polort $POLORT                                                             \
				-num_stimts 0                                                         \
				-fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
				-x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
				-fitts fitts.$subjstr                                                 \
				-errts errts.${subjstr}                                               \
				-x1D_stop                                                             \
				-bucket stats.$subjstr

			# -- use 3dTproject to project out regression matrix --
			#    (make errts like 3dDeconvolve, but more quickly)
			3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
					   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
					   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.tproject

			# if 3dDeconvolve fails, terminate the script
			if ( $status != 0 ) then
				echo '---------------------------------------'
				echo '** 3dDeconvolve error, failing...'
				echo '   (consider the file 3dDeconvolve.err)'
				exit
			endif

			# -- use 3dTproject to project out regression matrix --
			#    (make errts like 3dDeconvolve, but more quickly)
			3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
					   -mask mask_epi_anat.$subj+tlrc                                 \
					   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
					   -dsort Local_FSWMe_rall+tlrc                                   \
					   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.fanaticor

			endif

			# ================================ convert =================================

			# define niifile
			set niifile = $subj"_desc-s"$blur"preproc_bold.nii.gz"

			# do the AFNI to .nii conversion
			3dAFNItoNIFTI -prefix $derivoutdir/$niifile errts.${subjstr}.fanaticor+tlrc

			# remove temporary files
			\rm -f rm.*

		end

    endif

	cd $outdir

    ###########################################################################
    ############################## RESTING-STATE ##############################
    ###########################################################################

    set task = rest
    set subj = "$subject"_task-"$task"
    set POLORT = 5
	
	set output_dir = "$subj".results
	cd $output_dir
	set runs = (`count -digits 2 1 2`)

	################################### NO BLUR	################################### 

	set blur = 0
	set subjstr = "${subj}"_s"$blur"

    if ($run_rest == true) then

		# ================================= scale ==================================
		# scale each voxel time series to have a mean of 100
		# (be sure no negatives creep in)
		# (subject to a range of [0,200])
		foreach run ( $runs )
			3dTstat -prefix rm.mean_r$run pb03.$subj.r$run.volreg+tlrc
			3dcalc -a pb03.$subj.r$run.volreg+tlrc -b rm.mean_r$run+tlrc \
			       -c mask_epi_anat.$subj+tlrc                           \
			       -expr 'c * min(200, a/b*100)*step(a)*step(b)'         \
			       -prefix pb00.$subjstr.r$run.scale
		end

		# ================================ regress =================================
		# ------------------------------
		# run the regression analysis
		3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                         \
			-mask mask_epi_anat.$subj+tlrc                                        \
			-censor censor_${subj}_combined_2.1D                                  \
			-ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
			-ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
			-ortvec mot_demean.r01.1D mot_demean_r01                              \
			-ortvec mot_demean.r02.1D mot_demean_r02                              \
			-ortvec mot_deriv.r01.1D mot_deriv_r01                                \
			-ortvec mot_deriv.r02.1D mot_deriv_r02                                \
			-polort $POLORT                                                             \
			-num_stimts 0                                                         \
			-fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
			-x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
			-fitts fitts.$subjstr                                                 \
			-errts errts.${subjstr}                                               \
			-x1D_stop                                                             \
			-bucket stats.$subjstr

		# -- use 3dTproject to project out regression matrix --
		#    (make errts like 3dDeconvolve, but more quickly)
		3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
				   -mask mask_epi_anat.$subj+tlrc                                 \
				   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
				   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.tproject

		# if 3dDeconvolve fails, terminate the script
		if ( $status != 0 ) then
			echo '---------------------------------------'
			echo '** 3dDeconvolve error, failing...'
			echo '   (consider the file 3dDeconvolve.err)'
			exit
		endif

		# -- use 3dTproject to project out regression matrix --
		#    (make errts like 3dDeconvolve, but more quickly)
		3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
				   -mask mask_epi_anat.$subj+tlrc                                 \
				   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
				   -dsort Local_FSWMe_rall+tlrc                                   \
				   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.fanaticor

		# ================================ convert =================================

		# define niifile
		set niifile = $subj"_desc-s"$blur"preproc_bold.nii.gz"

		# do the AFNI to .nii conversion
		3dAFNItoNIFTI -prefix $derivoutdir/$niifile errts.${subjstr}.fanaticor+tlrc

		# remove temporary files
		\rm -f rm.*

		###################################   BLUR	################################### 

		foreach blur ( $kernels ) 

			set subjstr = "${subj}"_s"$blur"

			# ================================== blur ==================================
			# blur each volume of each run
			foreach run ( $runs )
				3dBlurToFWHM -FWHM $blur -mask mask_epi_anat.$subj+tlrc \
						     -input pb03.$subj.r$run.volreg+tlrc    \
						     -prefix pb00.$subjstr.r$run.blur 
			end

			# ================================= scale ==================================
			# scale each voxel time series to have a mean of 100
			# (be sure no negatives creep in)
			# (subject to a range of [0,200])
			foreach run ( $runs )
				3dTstat -prefix rm.mean_r$run pb00.$subjstr.r$run.blur+tlrc
				3dcalc -a pb00.$subjstr.r$run.blur+tlrc -b rm.mean_r$run+tlrc \
					   -expr 'min(200, a/b*100)*step(a)*step(b)'           \
					   -prefix pb00.$subjstr.r$run.scale
			end

			# ================================ regress =================================
			# ------------------------------
			# run the regression analysis
			3dDeconvolve -input pb00.$subjstr.r*.scale+tlrc.HEAD                         \
				-mask mask_epi_anat.$subj+tlrc                                        \
				-censor censor_${subj}_combined_2.1D                                  \
				-ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
				-ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
				-ortvec mot_demean.r01.1D mot_demean_r01                              \
				-ortvec mot_demean.r02.1D mot_demean_r02                              \
				-ortvec mot_deriv.r01.1D mot_deriv_r01                                \
				-ortvec mot_deriv.r02.1D mot_deriv_r02                                \
				-polort $POLORT                                                             \
				-num_stimts 0                                                         \
				-fout -tout -x1D X.s${blur}.xmat.1D -xjpeg X.s${blur}.jpg             \
				-x1D_uncensored X.s${blur}.nocensor.xmat.1D							  \
				-fitts fitts.$subjstr                                                 \
				-errts errts.${subjstr}                                               \
				-x1D_stop                                                             \
				-bucket stats.$subjstr

			# -- use 3dTproject to project out regression matrix --
			#    (make errts like 3dDeconvolve, but more quickly)
			3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
					   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
					   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.tproject

			# if 3dDeconvolve fails, terminate the script
			if ( $status != 0 ) then
				echo '---------------------------------------'
				echo '** 3dDeconvolve error, failing...'
				echo '   (consider the file 3dDeconvolve.err)'
				exit
			endif

			# -- use 3dTproject to project out regression matrix --
			#    (make errts like 3dDeconvolve, but more quickly)
			3dTproject -polort 0 -input pb00.$subjstr.r*.scale+tlrc.HEAD                 \
					   -mask mask_epi_anat.$subj+tlrc                                 \
					   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
					   -dsort Local_FSWMe_rall+tlrc                                   \
					   -ort X.s${blur}.nocensor.xmat.1D -prefix errts.${subjstr}.fanaticor

			endif

			# ================================ convert =================================

			# define niifile
			set niifile = $subj"_desc-s"$blur"preproc_bold.nii.gz"

			# do the AFNI to .nii conversion
			3dAFNItoNIFTI -prefix $derivoutdir/$niifile errts.${subjstr}.fanaticor+tlrc

			# remove temporary files
			\rm -f rm.*

		end

	endif

cd $outdir
							
end
