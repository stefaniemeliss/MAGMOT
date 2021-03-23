#!/bin/bash
source ~/.bashrc

module load afni19.3.03

path="/storage/shared/research/cinn/2018/MAGMOT"


# change directory to BIDS folder
BIDS_dir="$path"/MAGMOT_BIDS
cd $BIDS_dir

# declare deriv_dir
deriv_dir="$path"/derivatives

# declare ISC dir
ISC_root="$deriv_dir"/ISC
mkdir $ISC_root


# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
subjects=(sub-control001 sub-control002 sub-control003) # script development

# get length of subjects
num_subjects=${#subjects[@]}

# define task
task=rest

for (( s=0; s<${num_subjects}; s++));
    do
    # define variable subj_id, change directory to where the bold data of subj_id is saved and define file for subj_id
    subj_id=${subjects[$s]}
    s1_dir="$deriv_dir"/"$subj_id"/func
    cd $s1_dir
    s1_file=($(ls "$subj_id"_task-"$task"_desc-concat_bold.nii.gz))
    s1_files=($(ls "$subj_id"_task-"$task"_run-*))
    num_s1_files=${#s1_files[@]} # define the length of s1_files

    
    for (( t=s+1; t<${num_subjects}; t++));
        do

        # define variable subj_corr, change directory to where the bold data of subj_corr is saved and define file for subj_corr
        subj_corr=${subjects[$t]}
        s2_dir="$deriv_dir"/"$subj_corr"/func
        cd $s2_dir
        s2_file=($(ls "$subj_corr"_task-"$task"_desc-concat_bold.nii.gz))
        s2_files=($(ls "$subj_corr"_task-"$task"_run-*))
    	num_s2_files=${#s2_files[@]} # define the length of s2_files

        # create folder and compute ISC for whole time courrse
        ISC_dir="$ISC_root"/"$subj_id""$subj_corr"/
        mkdir $ISC_dir
        cd $ISC_dir

		# compute ISC for SME outcomes
       	for (( i=0; i<${num_s2_files}; i++)); do

		    # define files that should be correlated with one another
           	s1_file=${s1_files[$i]}
           	s2_file=${s2_files[$i]}

           	# define file prefix for ISC file
           	searchstring="${subj_corr}"
           	replacestring="${subj_id}""${subj_corr}"
           	ISC_prefix="ISC_""${s2_file/$searchstring/$replacestring}"

           	searchstring="_desc-preproc_bold.nii.gz"
           	replacestring="_z.nii.gz"
           	fisher_prefix="${ISC_prefix/$searchstring/$replacestring}"

            # compute pairwise ISC map
			if [ ! -f "$fisher_prefix" ]; then

	           	echo s1_file $s1_file
	           	echo s2_file $s2_file
           		echo ISC_SME_prefix $fisher_prefix
                echo ""
				# calculate ISC using pearson, do not detrend the data, save it in .nii, do Fisher-z transformation
				3dTcorrelate -pearson -polort -1 -Fisher -prefix $fisher_prefix  $s1_dir$s1_file $s2_dir$s2_file
			fi

		done

	done

done
