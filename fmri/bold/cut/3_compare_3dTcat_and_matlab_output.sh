#!/usr/bin/bash

# this code compares the files created with the Matlab and python code

# import modules
module load afni19.3.03

# define directories
path="/storage/shared/research/cinn/2018/MAGMOT"
cut_root="$path""/derivatives"
concat_root="$path""/derivatives/magictrickwatching/concat"

# define list of subjects
cd "$path"/MAGMOT_BIDS/ # change directory to BIDS folder

subjects=($(ls -d sub*)) # define subjects based on folder names in the BIDS directory
subjects=(sub-experimental005) # script development

num_subjects=${#subjects[@]} # get length of subject

# change directory
cd "$path""/derivatives/magictrickwatching/"
# create output file
printf "subject\tdiff_1_min\tdiff_1_max\tdiff_2_min\tdiff_2_max\tdiff_3_min\tdiff_3_max" > "$path""/derivatives/magictrickwatching/"matlab_3dTcat_diff.tsv #file header

for (( s=0; s<${num_subjects}; s++));
    do
    # define variable subj_id
    subj_id=${subjects[$s]}

	# define matlab file
	matlab_1="$concat_root"/"$subj_id"/"$subj_id"_task-magictrickwatching_run-1_bold_concat.nii.gz
	matlab_2="$concat_root"/"$subj_id"/"$subj_id"_task-magictrickwatching_run-2_bold_concat.nii.gz
	matlab_3="$concat_root"/"$subj_id"/"$subj_id"_task-magictrickwatching_run-3_bold_concat.nii.gz

	# define python file
	python_1="$cut_root"/"$subj_id"/func/"$subj_id"_task-magictrickwatching_run-1_desc-cut_bold.nii.gz
	python_2="$cut_root"/"$subj_id"/func/"$subj_id"_task-magictrickwatching_run-2_desc-cut_bold.nii.gz	
	python_3="$cut_root"/"$subj_id"/func/"$subj_id"_task-magictrickwatching_run-3_desc-cut_bold.nii.gz

	# calculate diff images
	3dcalc -a $matlab_1 -b $python_1 -expr 'a-b' -prefix 'diff_1.nii.gz'
	3dcalc -a $matlab_2 -b $python_2 -expr 'a-b' -prefix 'diff_2.nii.gz'
	3dcalc -a $matlab_3 -b $python_3 -expr 'a-b' -prefix 'diff_3.nii.gz'

	# determine min and max valkue in difference image
	min_1=$(3dBrickStat -min 'diff_1.nii.gz')
	max_1=$(3dBrickStat -max 'diff_1.nii.gz')
	min_2=$(3dBrickStat -min 'diff_2.nii.gz')
	max_2=$(3dBrickStat -max 'diff_2.nii.gz')
	min_3=$(3dBrickStat -min 'diff_3.nii.gz')
	max_3=$(3dBrickStat -max 'diff_3.nii.gz')

	# print values into file
	printf "\n$subj_id\t$min_1\t$max_1\t$min_2\t$max_2\t$min_3\t$max_3" >> "$path""/derivatives/magictrickwatching/"matlab_3dTcat_diff.tsv

	# delete all diff files
	rm diff*
 	
	
done
