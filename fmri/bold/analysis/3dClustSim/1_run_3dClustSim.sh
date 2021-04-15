#!/bin/bash
source ~/.bashrc

module load afni19.3.03

###############################################################################
########################## set up and run 3dClustSim ##########################
###############################################################################

# define path
DIR="/storage/shared/research/cinn/2018/MAGMOT"

# define derivatves dir
deriv_dir="$DIR"/derivatives

# set working directory
anal_root="$deriv_dir"/analysis

# define ROI dir
ROI_dir="$deriv_dir"/ROI_masks/output

# specify mask (created in previous script)
epi_mask=sample_label-dilatedGM_mask.nii.gz
gm_mask=sample_label-gm_mask.nii.gz

# define task
tasks=(magictrickwatching rest)
tasks=(rest)

# determine kernels
kernels=(s0 s4 s6 s8)
#kernels=(s0)

for task in "${tasks[@]}"; do

	task_root="$anal_root"/"$task"

	out_dir="$task_root"/3dClustSim
	cd $out_dir

	for blur in "${kernels[@]}"; do

		echo "$task $blur"

		# select first three columns of FWHMx output
		1d_tool.py -infile out_FWHMx_"$blur".1D -select_cols '0..2' -write out.1D

		# compute mean parameters based on subject-wise output from FWHMx
		new_params=$(3dTstat -mean -prefix stdout: out.1D\')

		rm out.1D # remove file

		echo $new_params

		# simulate cluster extent threshold
		#3dClustSim -acf $params -both -prefix ClustSim_"$task"_"$blur" -nxyz 64 76 64 -dxyz 3 3 3 -both -mask $ROI_dir/$epi_mask
		3dClustSim -acf $params -both -prefix ClustSim_"$task"_"$blur" -nxyz 64 76 64 -dxyz 3 3 3 -mask $ROI_dir/$epi_mask

	done

done
