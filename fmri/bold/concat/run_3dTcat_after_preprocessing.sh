#!/bin/bash

################################################################################
# concatenating task BOLD series after pre-processing
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
subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
subjects=(sub-control003)
#subjects=(sub-experimental016)

# define TSV file to read in
input=$DIR/"scripts/fmri/bold/concat/MAGMOT_inputForConcatenation.tsv"

# define search and replace string for file prefix
searchstring="_desc-preproc_bold.nii.gz"
replacestring="_desc-concat_bold.nii.gz"


#### FOR TESTING #####
#searchstring="_afniproc.nii.gz"
#replacestring="_concat.nii.gz"



# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

	# define output_dir
	out_dir=$deriv_dir/$subject/func/

	# load in pre-processed task file
	cd $out_dir
	file="$subject"_task-magictrickwatching_desc-preproc_bold.nii.gz

	# set variable to track appearance of first magic trick
	i=0
	first=1

	#### FOR TESTING #####
	#in_dir=$deriv_dir/magictrickwatching/concat/$subject/
	#cd $in_dir
	#file="$subject"_task-magictrickwatching_afniproc.nii.gz

	# read in file with information about scan duration
	{
	IGNORE_FL=T
	while read ID stim_file start_vol end_vol
		do

		# if the line is header set IGNORE_FL to F and skip line
		if [[ ${ID} = "ID" ]]; then
			IGNORE_FL=F
				continue
		fi

		# Ignore all lines until actual columns are found
		if [[ ${IGNORE_FL} = T  ]]; then
			continue
		fi

		# when the right information are available for the subject
		if [[ "$ID" == *"$subject"* ]]; then

			# update variable that tracks appearance of magictricks
			((i=i+1))

			# when we're dealing with the first magic trick
			if [[ $i -eq $first ]]; then

				echo "processing concatenation"

				# determine start and end vol for the first magic trick and save them in select
				select=`echo $start_vol..$end_vol`

			else

				# determine start and end vol for the next magic trick and save them in select_temp
				select_temp=`echo $start_vol..$end_vol`
				#update select to include start and end of all magic tricks processed so far
				select=`echo $select,$select_temp`

			fi


		fi

		done < "$input"
		}

	# once input file is processed
	echo "selected vols: $select"

	# 3dTcat to concatenate data using selected volumes
	prefix="${file/$searchstring/$replacestring}"
	3dTcat $file[$select] -dry -prefix $prefix -session $out_dir

done
