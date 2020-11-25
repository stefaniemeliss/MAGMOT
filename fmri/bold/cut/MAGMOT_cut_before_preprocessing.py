#!/usr/bin/python

# Cutting fMRI data
# This code is written to cut fMRI data along the 4th dimension, i.e. to reduce the number of volumes.
# This is necessary because the scanner did not stop automatically at the end of each block of the experiment, but had to be stopped manually.

# import necessary libraries
import pandas as pd
import nibabel as nib
import os
import re

# Load in information about scan duration - need to be in the same directory
# Please note that this file was downloaded using a seperate R script "1_download_from_OSF.R"
duration_file_path = os.getcwd()
duration_file_name = os.path.join(duration_file_path, 'MAGMOT_informationAboutScanDuration.tsv')
duration = pd.read_csv(duration_file_name, sep='\t')
duration.head()
#os.remove(duration_file_name) # delete file

# Create a list containing the subject names to loop through

# define path where the data set in BIDS format is stored
main_path = '/Users/stefaniemeliss/Dropbox/Reading/PhD/Magictricks/fmri_study' # local dir
main_path = '/storage/shared/research/cinn/2018/MAGMOT/' # remote dir
BIDS_path = os.path.join(main_path, 'MAGMOT_BIDS') # this has to be adapted as necessary
# list all folders in the BIDS directory
dir_list = os.listdir(BIDS_path)
# only use those that contain the string "sub"
subjects = [x for x in dir_list if "sub" in x]
subjects.sort() # sort subjects ascending
print(subjects)
print(len(subjects) == 50) # verify that the subject list has 50 elements
# subjects = ['sub-control001', 'sub-control002'] #debug

# define derivatives folder
derivatives_root = os.path.join(main_path, 'derivatives')
# create directory if it does not exist yet
if os.path.exists(derivatives_root) is not True:
    os.mkdir(derivatives_root)

### Use list to iterate through subjects and cut their fMRI data

for subj in subjects:

	# create folder for subject within derivatives
	subj_root = os.path.join(derivatives_root, subj)
	# create directory if it does not exist yet
	if os.path.exists(subj_root) is not True:
		os.mkdir(subj_root)
	# create func folder within subject within derivatives
	subj_func = os.path.join(subj_root, 'func')
	# create directory if it does not exist yet
	if os.path.exists(subj_func) is not True:
		os.mkdir(subj_func)

    # go into folder that contains the fMRI data
	data_path =  os.path.join(BIDS_path, subj, 'func')

    # create list that contains all task fMRI files
    # list all folders in the BIDS directory
	file_list = os.listdir(data_path)
    # only use those that contain the string "sub"
	task_files = [x for x in file_list if "magictrickwatching" in x]
	task_nii_files = [x for x in task_files if "nii" in x]
	task_nii_files.sort() # sort subjects ascending

    # define output path for each subject
	cut_path = os.path.join(subj_func, 'cut') # define output folder on subject level
	if os.path.exists(cut_path) is not True:
		os.mkdir(cut_path) # create outputfolder

    ####### create the next loop to iterate through through task_nii_files #######
	for scan_name in task_nii_files:
		print(scan_name)

        # define the file name of the nii file
		nii_filename = os.path.join(data_path, scan_name)

        # load in the nii file as img
		img = nib.load(nii_filename)

        # use scan_name to create scan_id: this will match the scan column in the duration dataframe
		scan_id = re.sub('_bold.nii.gz', '', scan_name)

        # get the duration of the scan_namen (i.e. scan_id) from the duration dataframe
        # find out which the corresponding row in the duration df is
		target_row = duration.loc[duration['scan'] == scan_id].index
		target_row_index = target_row.tolist() # transform Int64Index to integer
		print(target_row_index)
        # use target_row_index to get the actual duration in TR corresponding to the current scan_name
		target_duration = duration.loc[target_row_index, 'duration_run_TR']
		target_duration = target_duration.tolist() # transform Int64Index to integer
		target_duration_index = int(target_duration[0]) # we don't need to substract 1 because the stop slice is exclusive in python
		print(target_duration_index)

        # Slice nii file using .slicer
		cut = img.slicer[..., 0:target_duration_index] # drop the last couple of volumes

        # Save concatenated file
		cut_filename = re.sub('_bold.nii.gz', '_space-orig_desc-cut_bold.nii.gz', scan_name)
		cut_file = os.path.join(cut_path, cut_filename)
		nib.save(cut, cut_file)

