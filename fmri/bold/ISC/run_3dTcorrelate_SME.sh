#!/bin/bash
source ~/.bashrc

module load afni19.3.03

path="/storage/shared/research/cinn/2018/MAGMOT"


# change directory to BIDS folder
BIDS_dir="$path"/MAGMOT_BIDS
cd $BIDS_dir

# declare deriv_dir
deriv_dir="$path"/derivatives

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
subjects=(sub-control001 sub-control002 sub-control003) # script development

# get length of subjects
num_subjects=${#subjects[@]}

# define task
task=magictrickwatching

# set variable to track whether ISC_SME should be created or not
SME=F


for (( s=0; s<${num_subjects}; s++));
    do
    # define variable subj_id, change directory to where the bold data of subj_id is saved and define file for subj_id
    subj_id=${subjects[$s]}
    s1_dir="$deriv_dir"/"$subj_id"/func
    cd $s1_dir
    s1_file=($(ls "$subj_id"_task-"$task"_desc-concat_bold.nii.gz))
    
    for (( t=s+1; t<${num_subjects}; t++));
        do

        # define variable subj_corr, change directory to where the bold data of subj_corr is saved and define file for subj_corr
        subj_corr=${subjects[$t]}
        s2_dir="$deriv_dir"/"$subj_corr"/func
        cd $s2_dir
        s2_file=($(ls "$subj_corr"_task-"$task"_desc-concat_bold.nii.gz))

        # create folder and compute ISC for whole time courrse
        ISC_root="$deriv_dir"/ISC
        mkdir $ISC_root

        ISC_dir="$ISC_root"/"$subj_id""$subj_corr"/
        ISC_prefix="ISC_""$subj_id""$subj_corr""_task-"$task"_z.nii.gz"

        mkdir $ISC_dir
        cd $ISC_dir

		if [ ! -f "$ISC_prefix" ]; then
        	echo s1_file $s1_file
        	echo s2_file $s2_file
			echo ISC_prefix $ISC_prefix

			# calculate ISC using pearson, do not detrend the data, save it in .nii, do Fisher-z transformation
        	3dTcorrelate -pearson -polort -1 -Fisher -prefix $ISC_prefix $s1_dir$s1_file $s2_dir$s2_file

		fi


        # compute ISC maps according the ISC-SME principle
        if [[ ${SME} = T ]]; then

            # go into the concat pair folder
            concat_pair_dir="$path"/derivatives/magictrickwatching/concat/"$subj_id""$subj_corr"/
            cd $concat_pair_dir
            SME_s1=($(ls *_"$subj_id"*_SME_afniproc_bold.nii.gz))
            SME_s1=($(ls *_"$subj_id"*_SME_afniproc_bold.nii.gz))
            # echo SME_s1 $SME_s1
            num_SME_s1=${#SME_s1[@]} # define the length of SME_s1
            # echo num_SME_s1 $num_SME_s1
            SME_s2=($(ls *_"$subj_corr"*_SME_afniproc_bold.nii.gz))
            SME_s2=($(ls *_"$subj_corr"*_SME_afniproc_bold.nii.gz))
            # echo SME_s2 $SME_s2
            num_SME_s2=${#SME_s2[@]} # define the length of SME_s2
            # echo num_SME_s2 $num_SME_s2

		    # compute ISC for SME outcomes
           	for (( i=0; i<${num_SME_s2}; i++)); do

		        # define files that should be correlated with one another
               	SME_s1_file=${SME_s1[$i]}
               	echo SME_s1_file $SME_s1_file
               	SME_s2_file=${SME_s2[$i]}
               	echo SME_s2_file $SME_s2_file

               	# define file prefix for ISC_SME file
               	searchstring="${subj_id}${subj_corr}_${subj_id}_"
               	replacestring="_"
               	ISC_SME_prefix="ISC_""${subj_id}""${subj_corr}""${SME_s1_file/$searchstring/$replacestring}"
               		
			    searchstring=".nii.gz"
               	replacestring="_z.nii.gz"
               	fisher_prefix="${ISC_SME_prefix/$searchstring/$replacestring}"

			    # go into ISC folder 
            	cd $ISC_dir

			    if [ ! -f "$fisher_prefix" ]; then

               		echo ISC_SME_prefix $fisher_prefix
				    # calculate ISC using pearson, do not detrend the data, save it in .nii, do Fisher-z transformation
				    3dTcorrelate -pearson -polort -1 -Fisher -prefix $fisher_prefix  $concat_pair_dir$SME_s1_file $concat_pair_dir$SME_s2_file
			    fi

		    done

		fi

	done
done
