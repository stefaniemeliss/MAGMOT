#!/bin/bash
source ~/.bashrc

# load modules
module load afni19.3.03
# load anaconda - this is necessary to have all python libraries available
module load anaconda3

#define path
path="/storage/shared/research/cinn/2018/MAGMOT"

# define directories
BIDS_dir="$path"/MAGMOT_BIDS
deriv_dir="$path"/derivatives
software_dir="$path"/software/pyfMRIqc-master

# change directory to BIDS folder
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))

#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006) # script development
#subjects=(sub-control037) # script development
#subjects=(sub-control001)

# create file to save censor values
#printf "subject\tscan\tvolumes\tthreshold_02mm\tfraction_02mm\tthreshold_03mm\tfraction_03mm\tthreshold_05mm\tfraction_05mm\tthreshold_1mm\tfraction_1mm" > "$path"/derivatives/pyfMRIqc/afni_volreg_quick_censor_count.tsv #file header
printf "subject\tscan\tvolumes\tthreshold_02\tfraction_02\tthreshold_03\tfraction_03\tthreshold_05\tfraction_05\tthreshold_1\tfraction_1" > "$path"/derivatives/pyfMRIqc/afni_volreg_censor_count.tsv #file header


# for each subject in the subjects array
for subject in "${subjects[@]}"; do
	echo $subject

	# define all necessary directories
	fildirTASK=$deriv_dir/$subject/func	
	fildirREST=$BIDS_dir/$subject/func
	fildirQC=$deriv_dir/pyfMRIqc/$subject	
	mkdir $fildirQC

	# copy task files
	cd $fildirTASK
	filesTASK=($(ls -d *magictrickwatching*cut_bold.nii.gz))
	for fileTASK in "${filesTASK[@]}"; do
		searchstring="nii.gz"
		replacestring="nii.gz"	
		fileTASK_new="${fileTASK/$searchstring/$replacestring}"	
		3dcopy $fildirTASK/$fileTASK $fildirQC/$fileTASK_new
	done	

	# copy rest files
	cd $fildirREST
	filesREST=($(ls -d *rest*.nii.gz))
	for fileREST in "${filesREST[@]}"; do
		searchstring="nii.gz"
		replacestring="nii.gz"	
		fileREST_new="${fileREST/$searchstring/$replacestring}"	
		3dcopy $fildirREST/$fileREST $fildirQC/$fileREST_new
	done

	# create files array with nifti files
	cd $fildirQC
	files=($(ls -d sub*.nii*))

	# loop through the files that have been copied over
	for file in "${files[@]}"; do

		cd $fildirQC

		#########################
		### MOTION CORRECTION ###
		#########################

		# define file prefix for motion corrected files
		searchstring="sub"
		replacestring_3d="volreg_sub"	
		replacestring_1d="motion_sub"
		replacestring_plot="plot_sub"
		searchending=".nii.gz"
		replaceending=""
		volreg_prefix="${file/$searchstring/$replacestring_3d}"	
		file_suffix="${file/$searchending/$replaceending}" 
		motionfile_prefix="${file_suffix/$searchstring/$replacestring_1d}".1D	
		plot_prefix="${file_suffix/$searchstring/$replacestring_plot}"

		# run motion correction
		3dvolreg -verbose -prefix $volreg_prefix -1Dfile $motionfile_prefix  $file 

		# look at the 1dplot
		1dplot -volreg -sepsc1 -png $plot_prefix $motionfile_prefix &

		# extract number of censored TRs for different thresholds
		# using -quick_censor_count is akin to -censor_motion, but it only outputs the number of TRs that would be censored
        # This option essentially replaces these: -derivative -demean -collapse_cols euclidean_norm -censor_prev_TR -verb 0 -show_censor_count -moderate_mask 0 LIMIT
		# but censor_prev_TR is not needed
		# thresh_02=$(1d_tool.py -infile $motionfile_prefix -quick_censor_count 0.2)
		# thresh_03=$(1d_tool.py -infile $motionfile_prefix -quick_censor_count 0.3)
		# thresh_05=$(1d_tool.py -infile $motionfile_prefix -quick_censor_count 0.5)
		# thresh_1=$(1d_tool.py -infile $motionfile_prefix -quick_censor_count 1)

		thresh_02=$(1d_tool.py -infile $motionfile_prefix -derivative -demean -collapse_cols euclidean_norm -verb 0 -show_censor_count -moderate_mask 0 0.2)
		thresh_03=$(1d_tool.py -infile $motionfile_prefix -derivative -demean -collapse_cols euclidean_norm -verb 0 -show_censor_count -moderate_mask 0 0.3)
		thresh_05=$(1d_tool.py -infile $motionfile_prefix -derivative -demean -collapse_cols euclidean_norm -verb 0 -show_censor_count -moderate_mask 0 0.5)
		thresh_1=$(1d_tool.py -infile $motionfile_prefix -derivative -demean -collapse_cols euclidean_norm -verb 0 -show_censor_count -moderate_mask 0 1)


		# calculate percentage of censored volumes
		num_vols=$(3dinfo -nt $file)
		frac_02=$(bc <<<"scale=5; $thresh_02 / $num_vols")
		frac_03=$(bc <<<"scale=5; $thresh_03 / $num_vols")
		frac_05=$(bc <<<"scale=5; $thresh_05 / $num_vols")
		frac_1=$(bc <<<"scale=5; $thresh_1 / $num_vols")


		# print values to file
		# printf "\n$subject\t$file\t$num_vols\t$thresh_02\t$frac_02\t$thresh_03\t$frac_03\t$thresh_05\t$frac_05\t$thresh_1\t$frac_1" >> "$path"/derivatives/pyfMRIqc/afni_volreg_quick_censor_count.tsv
		printf "\n$subject\t$file\t$num_vols\t$thresh_02\t$frac_02\t$thresh_03\t$frac_03\t$thresh_05\t$frac_05\t$thresh_1\t$frac_1" >> "$path"/derivatives/pyfMRIqc/afni_volreg_censor_count.tsv

		##########################
		###### CREATE MASK  ######
		##########################

		# define file prefix for mask
		searchstring="sub"
		replacestring_mask="mask_sub"	
		mask_prefix="${file/$searchstring/$replacestring_mask}"	

		# use 3dAutomask to perform skull-stripping and masking on EPI data (default in pre-processing)
		3dAutomask -prefix $mask_prefix $file

		##########################
		###### RUN pyfMRIqc ######
		##########################

		# change directory to pyfMRIqc-master
		#cd /storage/shared/research/cinn/2018/MAGMOT/software/pyfMRIqc-fix-rm
		#cd /storage/shared/research/cinn/2018/MAGMOT/software/pyfMRIqc-master

		# define input parameters
		func_nift_file=$fildirQC/"${file}"
		SNR_voxel_perc=25
		#mask_threshold=200
		motion_file=$fildirQC/"${motionfile_prefix}"
		mask_nift_file=$fildirQC/"${mask_prefix}"

		#python pyfMRIqc.py -n $func_nift_file -s $SNR_voxel_perc -m $motion_file -k $mask_nift_file
		python $software_dir/pyfMRIqc.py -n $func_nift_file -s $SNR_voxel_perc -m $motion_file -k $mask_nift_file

	done 

	# remove nifti files
	rm *.nii.gz
	
done






