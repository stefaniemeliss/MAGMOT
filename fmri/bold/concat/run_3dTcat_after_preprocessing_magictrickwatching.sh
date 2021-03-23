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
#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
#subjects=(sub-control001)
#subjects=(sub-experimental016 sub-control045)

# define TSV file to read in
input=$DIR/"scripts/fmri/bold/concat/MAGMOT_inputForConcatenation.tsv"

# define search and replace string for file prefix
searchstring="preproc_bold.nii.gz"
replacestring="concat_bold.nii.gz"

searchstring="nosmooth_bold.nii.gz"
replacestring="nosmoothconcat_bold.nii.gz"

searchstring="smoothed_bold.nii.gz"
replacestring="smoothedconcat_bold.nii.gz"

#searchstring="smoothed_bold.nii.gz"
#replacestring="smoothedconcat_bold.nii.gz"

# define task
task=magictrickwatching

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
	niifile="$subject"_task-"$task"_desc-"$searchstring"

	# do the AFNI to .nii conversion
	3dAFNItoNIFTI -prefix $outdir/$niifile $file

    # remove file to save disc sapce
    rm $file_string


    ############### concatenate pre-processed file ###############

	# load in pre-processed task file
	cd $outdir

	# set variable to track appearance of first magic trick
	i=0
	first=1

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
	prefix="${niifile/$searchstring/$replacestring}"
	3dTcat $niifile[$select] -prefix $prefix -session $outdir


    # extract final file length
    3dinfo $prefix > info.txt
    numvol=$(awk 'NR==18 {print $6}' info.txt)
    echo "final number of vols: $numvol"
    rm info.txt

    # delete file to save disk space
    rm $niifile

done
