#!/bin/bash

################################################################################
# copy BIDS data and cut task BOLD series before pre-processing
################################################################################

# define path
DIR="/storage/shared/research/cinn/2018/MAGMOT"

# define rawdata dir
out_dir="$DIR"/rawdata
mkdir $out_dir

# change directory to BIDS folder
BIDS_dir="$DIR"/MAGMOT_BIDS
# change directory to the raw NIFI files
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006)
#subjects=(sub-control003)
#subjects=(sub-experimental016)

# define TSV file to read in
input=$DIR/"scripts/fmri/bold/cut/MAGMOT_inputForCutting.tsv"

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

	# define directories
	BIDS_dir_subj=$BIDS_dir/$subject

	# create output directories
	out_dir_subj=$out_dir/$subject
	mkdir $out_dir_subj
	mkdir $out_dir_subj/anat
	mkdir $out_dir_subj/fmap
	mkdir $out_dir_subj/func

	# copy tsv file
	cp $BIDS_dir_subj/*.tsv $out_dir_subj/

 	# copy fmap and anat folder
	cp -R $BIDS_dir_subj/anat $out_dir_subj
	cp -R $BIDS_dir_subj/fmap $out_dir_subj

	# copy rest EPI files
	cp $BIDS_dir_subj/func/*rest* $out_dir_subj/func/

	# copy json EPI files
	cp $BIDS_dir_subj/func/*.json $out_dir_subj/func/

	# list task files
	func_dir=$BIDS_dir_subj/func
	cd $func_dir
	files=($(ls  *magictrickwatching*nii.gz))

    # cut task files
	for file in "${files[@]}"; do

			# read in file with information about scan duration
			{
			IGNORE_FL=T
			while read ID scan duration_run_seconds duration_scan_seconds duration_run_TR duration_scan_TR
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
				if [[ "$file" == *"$scan"* ]]; then
					echo "It's there in $file"
					echo "scan $scan should have $duration_run_TR volumes"

					# determine the first volume 3dTcat should concatenate (this is constant)
					start_vol=0 #note that afni starts at 0
					end_vol=`expr $duration_run_TR - 1` # subtracting 1 because afni starts at 0
					echo "first vol: $start_vol and last vol: $end_vol"

					# 3dTcat to concatenate data using start_vol & end_vol to determine first and last volume to select
					3dTcat $file[$start_vol..$end_vol] -prefix $file -session $out_dir_subj/func
				fi

			done < "$input"
			}

	done
done

# at the end copy data set and participant files
cp $BIDS_dir/* $out_dir/

