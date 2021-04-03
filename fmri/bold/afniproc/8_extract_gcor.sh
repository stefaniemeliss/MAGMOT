#!/bin/bash

################################################################################
# extract global correlation
################################################################################


# define path
DIR="/storage/shared/research/cinn/2018/MAGMOT"

# change directory to BIDS folder
BIDS_dir="$DIR"/rawdata
# change directory to the raw NIFI files
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))

# sort array
subjects=($(echo ${subjects[*]}| tr " " "\n" | sort -n))

#subjects=(sub-control001)

# define derivatives directory
deriv_dir="$DIR"/derivatives

# create header of output file
printf "File\tgcor_magictrickwatching" > "$deriv_dir"/afniproc/gcor_magictrickwatching.txt #file header
printf "File\tgcor_rest" > "$deriv_dir"/afniproc/gcor_rest.txt #file header

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	# create output folder
	out_dir=$deriv_dir/$subject/func
    cd $out_dir

    ############## MAGICTRCIKWATCHING ##############

    task=magictrickwatching

    errts="$subject"_task-"$task"_desc-fullpreproc_bold.nii.gz
    mask="$subject"_task-"$task"_space-MNI152NLin2009cAsym_label-dilatedGM_mask.nii.gz

    # ---------------------------------------------------
    # compute and store GCOR (global correlation average)
    # (sum of squares of global mean of unit errts)
    3dTnorm -norm2 -prefix rm.errts.unit $errts
    3dmaskave -quiet -mask $mask rm.errts.unit+tlrc            \
              > mean.errts.unit.1D
    3dTstat -sos -prefix - mean.errts.unit.1D\' > out.gcor.1D
    echo "-- GCOR = `cat out.gcor.1D`"

    gcor=$(awk '{print $0}' out.gcor.1D)

    rm out.gcor.1D
    rm *errts.unit*

    printf "\n$errts\t$gcor" >> "$deriv_dir"/afniproc/gcor_magictrickwatching.txt


    ############## RESTING-STATE ##############

    task=rest

    errts="$subject"_task-"$task"_desc-fullpreproc_bold.nii.gz
    mask="$subject"_task-"$task"_space-MNI152NLin2009cAsym_label-dilatedGM_mask.nii.gz

    # ---------------------------------------------------
    # compute and store GCOR (global correlation average)
    # (sum of squares of global mean of unit errts)
    3dTnorm -norm2 -prefix rm.errts.unit $errts
    3dmaskave -quiet -mask $mask rm.errts.unit+tlrc            \
              > mean.errts.unit.1D
    3dTstat -sos -prefix - mean.errts.unit.1D\' > out.gcor.1D
    echo "-- GCOR = `cat out.gcor.1D`"

    gcor=$(awk '{print $0}' out.gcor.1D)

    rm out.gcor.1D
    rm *errts.unit*

    printf "\n$errts\t$gcor" >> "$deriv_dir"/afniproc/gcor_rest.txt


done


