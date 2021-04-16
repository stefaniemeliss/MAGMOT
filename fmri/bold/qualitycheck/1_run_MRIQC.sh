#!/bin/bash
source ~/.bashrc

#define path
path="/Users/stefaniemeliss/Dropbox/Reading/PhD/Magictricks/fmri_study"

# change directory to BIDS folder
cd "$path"/rawdata/

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))

#subjects=(sub-control001 sub-control002 sub-control003 sub-experimental004 sub-experimental005 sub-experimental006) # script development
#subjects=(sub-control041) # script development

# for each subject in the subjects array
for subject in "${subjects[@]}"; do
	echo $subject

	# run mriqc in docker container
	# on subject level
	# verbose reports to increase number of images in the html report
	# threshold for frame wise displacement 0.5
	# output workflow graph
	# perform despiking during motion correction
	# perform slice timing correction

  docker run -it --rm -v "$path"/rawdata:/data:ro -v "$path"/derivatives/mriqc:/out poldracklab/mriqc:latest /data /out participant --participant-label $subject --verbose-reports --write-graph --ica --fd_thres 0.5 --despike --correct-slice-timing

done

# create group report
docker run -it --rm -v "$path"/rawdata:/data:ro -v "$path"/derivatives/mriqc:/out poldracklab/mriqc:latest /data /out group

# upload report to osfclient
pip install osfclient
projectid="fhqb7"
osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force ~/Dropbox/Reading/PhD/Magictricks/fmri_study/derivatives/mriqc/group_bold.tsv osfstorage/data/fmri/group_bold.tsv
osf -p $projectid -u stefanie.meliss@pgr.reading.ac.uk upload --force ~/Dropbox/Reading/PhD/Magictricks/fmri_study/derivatives/mriqc/group_T1w.tsv osfstorage/data/fmri/group_T1w.tsv
