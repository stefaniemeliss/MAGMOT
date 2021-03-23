#!/bin/bash
source ~/.bashrc

module load afni19.3.03

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"

# define task
task=rest

# define directories
BIDS_dir="$path"/MAGMOT_BIDS
deriv_dir=$path/derivatives
ROI_dir="$deriv_dir"/ROI_masks/output
anal_dir=$deriv_dir/analysis/"$task"
CS_dir=$anal_dir/3dClustSim
FC_dir=$anal_dir/FC
code_dir=$path/scripts/fmri/bold/analysis/FC

# change directory to BIDS folder
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001)
# sort array
subjects=($(echo ${subjects[*]}| tr " " "\n" | sort -n))

preproc=(smoothed nosmooth)
preproc=(preproc nosmooth)

for smooth in "${preproc[@]}"; do

    # define ROI masks
    VTA_fHPC=$FC_dir/VTA_RSFC_"$smooth"_run-1_pearson_rewardHPC.nii.gz
    VTA_aHPC=$FC_dir/VTA_RSFC_"$smooth"_run-1_pearson_aHPC.nii.gz

VTA_RSFC_nosmooth_run-1_pearson_aHPC.nii.gz

    fHPC=$ROI_dir/MNI_res-epi_label-HPC_desc-rewardsensitiv_mask.nii.gz
    aHPC=$ROI_dir/MNI_res-epi_label-aHPC_mask.nii.gz

    # create file to save correlation values
    out_FC_fHPC=FC_"$smooth"_VTA-fHPC.txt
    out_FC_aHPC=FC_"$smooth"_VTA-aHPC.txt

    # open file for correlation output
    printf "BIDS\tRSFC_"$smooth"_fHPC_run-1\tRSFC_"$smooth"_fHPC_run-2\ttaskFC_"$smooth"_fHPC" > "$FC_dir"/$out_FC_fHPC #file header
    printf "BIDS\tRSFC_"$smooth"_aHPC_run-1\tRSFC_"$smooth"_aHPC_run-2\ttaskFC_"$smooth"_aHPC" > "$FC_dir"/$out_FC_aHPC #file header

    # define strings for file prefix
    replacestring="maskave_"

    # for each subject in the subjects array
    for subject in "${subjects[@]}"; do

	    echo $subject

	    # go to directory with concat files
	    concat_dir=$deriv_dir/$subject/func
	    cd $concat_dir

	    # PRE-PROCESSED AND CONCATENATED RESTING FILES
        if [[ "$smooth" == *"nosmooth"* ]]; then
            inputs=("$subject"_task-rest_run-1_desc-"$smooth"_bold.nii.gz "$subject"_task-rest_run-2_desc-"$smooth"_bold.nii.gz "$subject"_task-magictrickwatching_desc-"$smooth"concat_bold.nii.gz)
        else
            inputs=("$subject"_task-rest_run-1_desc-"$smooth"_bold.nii.gz "$subject"_task-rest_run-2_desc-"$smooth"_bold.nii.gz "$subject"_task-magictrickwatching_desc-concat_bold.nii.gz)
        fi

	    # loop over the preprocessed resting state files and extract average time course
	    for (( f=0; f<${#inputs[@]}; f++)); do

		    # define file_id
	        file=${inputs[$f]}

            echo $file

            # define search string based on whether its rest or task data
            if [[ "$file" == *"rest"* ]]; then
                searchstring="desc-"$smooth"_bold.nii.gz"
            else
                if [[ "$smooth" == *"nosmooth"* ]]; then
                    searchstring="desc-"$smooth"concat_bold.nii.gz"
                else
                    searchstring="desc-concat_bold.nii.gz"
                fi
            fi

            echo $searchstring

		    # define prefix
		    maskave_prefix="${file/$searchstring/$replacestring}"
		    VTA_fHPC_prefix="$maskave_prefix"VTA_fHPC.txt
		    VTA_aHPC_prefix="$maskave_prefix"VTA_aHPC.txt
		    fHPC_prefix="$maskave_prefix"fHPC.txt
		    aHPC_prefix="$maskave_prefix"aHPC.txt
		    FC_fHPC_prefix="$maskave_prefix"FC_fHPC.txt
		    FC_aHPC_prefix="$maskave_prefix"FC_aHPC.txt

		    # extract ROI average time course: for each volume, the average of voxel is computed in saved in txt file
		    3dmaskave -quiet -mask $VTA_fHPC $file > $VTA_fHPC_prefix # time course of VTA (determined with reward-sensitive HPC)
		    3dmaskave -quiet -mask $VTA_aHPC $file > $VTA_aHPC_prefix # time course of VTA (determined with anterior HPC)
		    3dmaskave -quiet -mask $fHPC $file > $fHPC_prefix # time course of reward-sensitive HPC
		    3dmaskave -quiet -mask $aHPC $file > $aHPC_prefix # time course of anterior HPC

		    # compute pearson correlation between VTA and fHPC
		    1dCorrelate -Pearson $VTA_fHPC_prefix $fHPC_prefix > $FC_fHPC_prefix
		    # grep cor value from file
		    pearson_f=$(awk 'NR==4 {print $3}' $FC_fHPC_prefix)

		    # compute pearson correlation between VTA and aHPC
		    1dCorrelate -Pearson $VTA_aHPC_prefix $aHPC_prefix > $FC_aHPC_prefix
		    # grep cor value from file
		    pearson_a=$(awk 'NR==4 {print $3}' $FC_aHPC_prefix)

            # add correlation to table
		    if [[ $f == 0 ]]; then
			    printf "\n$subject\t$pearson_f" >> "$FC_dir"/$out_FC_fHPC # cor for run-1
			    printf "\n$subject\t$pearson_a" >> "$FC_dir"/$out_FC_aHPC # cor for run-1
		    else
			    printf "\t$pearson_f" >> "$FC_dir"/$out_FC_fHPC # cor for run-2 + task
			    printf "\t$pearson_a" >> "$FC_dir"/$out_FC_aHPC # cor for run-2 + task
		    fi

	    done

    done

    # remove all the files that were created to calculate RSFC
    #rm $deriv_dir/*/func/*maskave_*

done


