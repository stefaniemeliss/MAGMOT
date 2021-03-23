#!/bin/bash
source ~/.bashrc

# load module
module load afni19.3.03
module load anaconda3


path="/storage/shared/research/cinn/2018/MAGMOT"

# define directories
deriv_dir="$path"/derivatives
afniproc_dir="$deriv_dir"/afniproc
 

# change directory to ROI folder
ROI_root="$deriv_dir"/ROI_masks
ROI_dir="$ROI_root"/output
mkdir $ROI_dir
input_dir="$ROI_root"/input
anat_dir="$ROI_root"/anatomical_masks

cd $ROI_dir

# define which template to use and where to find them
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`

# copy MNI template that has been used to normalise data to ROI_dir
3dcopy $template_path/$template $ROI_dir/$template

###### create average EPI masks ######

# define prefix
epi_ave_rest=sample_task-rest_label-automask_average.nii.gz
epi_mask_rest=sample_task-rest_label-automask_mask.nii.gz

epi_ave_mtw=sample_task-magictrickwatching_label-automask_average.nii.gz
epi_mask_mtw=sample_task-magictrickwatching_label-automask_mask.nii.gz

# create average EPI masks restingstate: Combine individual EPI automasks into a group mask
mask_ea=mask_epi_anat_010.*+tlrc.HEAD

3dMean -datum float -prefix $ROI_dir/$epi_ave_rest  "$afniproc_dir"/sub*/*_task-rest.results/$mask_ea
3dcalc -datum byte -prefix $ROI_dir/$epi_mask_rest -a $ROI_dir/$epi_ave_rest -expr 'step(a-0.499)'

mask_ea=mask_epi_anat.*+tlrc.HEAD

# create average EPI masks magic trick watching: Combine individual EPI automasks into a group mask
3dMean -datum float -prefix $ROI_dir/$epi_ave_mtw "$afniproc_dir"/sub*/*_task-magictrickwatching.results/$mask_ea
3dcalc -datum byte -prefix $ROI_dir/$epi_mask_mtw -a $ROI_dir/$epi_ave_mtw -expr 'step(a-0.499)'

###### create GM MASK ######

# use reference template volume 5 for @SSwarper https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/template_atlas/sswarper_base.html

#[0] = skull-stripped template brain volume
#[1] = skull-on template brain volume
#[2] = weight mask for nonlinear registration, with the brain given greater weight than the skull
#[3] = binary mask for the brain
#[4] = binary mask for gray matter plus some CSF (slightly dilated)
#      -- this volume is not used in this script
#      -- it is intended for use in restricting FMRI analyses to the 'interesting' parts of the brain
#      -- this mask should be resampled to your EPI spatial resolution (see program 3dfractionize), and then combined with a mask from your experiment reflecting your EPI brain coverage (see program 3dmask_tool).

# define prefix
gm_t1=MNI_res-anat_label-GM_mask.nii.gz
gm_epi=MNI_res-epi_label-GM_mask.nii.gz

# select 5th volume
3dTcat $template[4] -prefix $gm_t1

# resample mask to EPI resolution
3dresample -master $epi_ave_mtw -input $gm_t1 -prefix $gm_epi

# define prefix
rest_union=sample_task-rest_label-GM_desc-MNIunion_mask.nii.gz
mtw_union=sample_task-magictrickwatching_label-GM_desc-MNIunion_mask.nii.gz

# combine GM mask with EPI brain coverage 
3dmask_tool -input $gm_epi $epi_mask_rest -prefix $rest_union -union
3dmask_tool -input $gm_epi $epi_mask_mtw -prefix $mtw_union -union

# define prefix
rest_inter=sample_task-rest_label-GM_desc-MNIinter_mask.nii.gz
mtw_inter=sample_task-magictrickwatching_label-GM_desc-MNIinter_mask.nii.gz

# combine GM mask with EPI brain coverage 
3dmask_tool -input $gm_epi $epi_mask_rest -prefix $rest_inter -inter
3dmask_tool -input $gm_epi $epi_mask_mtw -prefix $mtw_inter -inter

