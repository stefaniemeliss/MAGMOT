#!/bin/bash

################################################################################
# concatenating rest BOLD series after pre-processing
################################################################################

# load in AFNI module
module load afni19.3.03

# define path
DIR="/storage/shared/research/cinn/2018/MAGMOT"

# define derivatves dir
deriv_dir="$DIR"/derivatives

# change directory to BIDS folder
BIDS_dir="$DIR"/MAGMOT_BIDS
# change directory to the raw NIFI files
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
#subjects=(sub-control001)
#subjects=(sub-experimental016)

# define search and replace string for file prefix
searchstring="desc-preproc_bold.nii.gz"
replacestring="desc-concat_bold.nii.gz"

replacestring1="run-1_desc-preproc_bold.nii.gz"
replacestring2="run-2_desc-preproc_bold.nii.gz"


searchstring="desc-nosmooth_bold.nii.gz"

replacestring1="run-1_desc-nosmooth_bold.nii.gz"
replacestring2="run-2_desc-nosmooth_bold.nii.gz"

searchstring="desc-smoothed_bold.nii.gz"

replacestring1="run-1_desc-smoothed_bold.nii.gz"
replacestring2="run-2_desc-smoothed_bold.nii.gz"

# define task
task=rest


# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

    ############### convert AFNI to NIFTI ###############

	# go to the $subject.magictrickwatching_perRun.results folder
	procdir=$deriv_dir/afniproc/$subject/"$subject"_task-"$task".results
	procdir=$deriv_dir/afniproc/$subject/"$subject"_task-"$task".results_old

	# cd into subject folder
	cd $procdir

	# define files: final preprocessing output is errts.$subjects.magictrickwatching_perRun.fanaticor
    file_string=errts."$subject"_task-"$task".fanaticor+tlrc*
    file_string=errts_010."$subject"_task-"$task".fanaticor+tlrc*
    file_string=errts_nosmooth."$subject"_task-"$task".fanaticor+tlrc*
    file_string=errts_smoothed."$subject"_task-"$task".fanaticor+tlrc*

	file=($(ls -d $file_string))
	#file=($(ls -d errts_smoothed."$subject"_task-"$task".fanaticor+tlrc*))

	# define concat directory (output for 3dcopy
	outdir=$deriv_dir/$subject/func

	# define niifile
	niifile="$subject"_task-"$task"_"$searchstring"

	# do the AFNI to .nii conversion
	3dAFNItoNIFTI -prefix $outdir/$niifile $file

    # remove file to save disc sapce
    #rm $file_string


    ############### concatenate pre-processed file ###############

	# load in pre-processed task file
	cd $outdir

	# define tcat prefix
	tcat_prefix1="${niifile/$searchstring/$replacestring1}"
	tcat_prefix2="${niifile/$searchstring/$replacestring2}"

	# create run-1 (first 300 vols)
	if [ ! -f "$tcat_prefix1" ]; then
		echo $tcat_prefix1
		3dTcat -prefix $tcat_prefix1 $niifile[0..299]
	fi

	# create run-2 (second 300 vols)
	if [ ! -f "$tcat_prefix2" ]; then
		echo $tcat_prefix2
		3dTcat -prefix $tcat_prefix2 $niifile[300..599]
	fi

    # remove nii file to save disc space
    #rm $niifile

done
