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

############### create average EPI and gray matters masks ###############

#cd "$deriv_dir" # change directory to where pre-processed files are

task=both

# define prefix
epi_ave=sample_task-"$task"_label-automask_average.nii.gz
epi_mask=sample_label-dilatedGM_mask.nii.gz

gm_ave=sample_task-"$task"_label-gm_average.nii.gz
gm_mask=sample_label-gm_mask.nii.gz

# define string for the masks that are created during pre-processing
mask_ea=*_space-MNI152NLin2009cAsym_label-dilatedGM_mask.nii.gz

# create average EPI masks restingstate: Combine individual EPI automasks into a group mask
3dMean -datum float -prefix $ROI_dir/$epi_ave  $deriv_dir/sub*/func/$mask_ea
3dcalc -datum float -prefix $ROI_dir/$epi_mask -a $ROI_dir/$epi_ave -expr 'ispositive(a-0.499)'

# remove average
rm $ROI_dir/$epi_ave

# define string for GM mask in EPI resolution and MNI space (note: this mask is the same for task and rest)
mask_gm=*_task-rest.results/follow_ROI_FSGMe+tlrc.BRIK.gz

# create average EPI masks restingstate: Combine individual EPI automasks into a group mask
3dMean -datum float -prefix $ROI_dir/$gm_ave  $deriv_dir/afniproc/sub-*/$mask_gm
3dcalc -datum float -prefix $ROI_dir/$gm_mask -a $ROI_dir/$gm_ave -expr 'ispositive(a-0.09)'

# remove average
rm $ROI_dir/$gm_ave


###### create average EPI masks ######

# define prefix
epi_ave_rest=sample_task-rest_label-automask_average.nii.gz
epi_mask_rest=sample_task-rest_label-automask_mask.nii.gz

epi_ave_mtw=sample_task-magictrickwatching_label-automask_average.nii.gz
epi_mask_mtw=sample_task-magictrickwatching_label-automask_mask.nii.gz

mask_ea=mask_epi_anat_010*+tlrc.HEAD

# create average EPI masks restingstate: Combine individual EPI automasks into a group mask
3dMean -datum float -prefix $ROI_dir/$epi_ave_rest  "$afniproc_dir"/sub*/*_task-rest.results/$mask_ea
3dcalc -datum byte -prefix $ROI_dir/$epi_mask_rest -a $ROI_dir/$epi_ave_rest -expr 'step(a-0.499)'

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

################################ create ROIs ################################

prefix=MNI_label-
suffix=_mask.nii.gz

### REWARD NETWORK ###

# note: this assumes that 1_run_atlaskit was executed already

# combine midbrain masks
3dcalc -a $anat_dir/"$prefix"VTA"$suffix" -b $anat_dir/"$prefix"SNc"$suffix" -c $anat_dir/"$prefix"SNr"$suffix" -expr 'a+b+c' -prefix $anat_dir/"$prefix"VTA_SN"$suffix"

# binarise SN/VTA mask
3dmask_tool -input $anat_dir/"$prefix"VTA_SN"$suffix" -prefix $anat_dir/"$prefix"VTASN"$suffix"

# remove SNc/SNr/VTA
rm $anat_dir/*SNc_*
rm $anat_dir/*SNr_*
rm $anat_dir/*VTA_*

### HIPPOCAMPUS ###

# use AFNI to extract HPC mask from Glasser et al 2016 atlas
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:L_Hippocampus -prefix $anat_dir/"$prefix"HPC_L"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:R_Hippocampus -prefix $anat_dir/"$prefix"HPC_R"$suffix"

# combine masks for both hemispheres and delete unilateral masks
3dcalc -a $anat_dir/"$prefix"HPC_L"$suffix" -b $anat_dir/"$prefix"HPC_R"$suffix" -expr 'a+b' -prefix $anat_dir/"$prefix"HPC"$suffix"
rm $anat_dir/*HPC_L* $anat_dir/*HPC_R*

# copy masks for anterior and posterior HPC to anatomical dir
3dcopy $input_dir/Gruber_HPC/0426rrewrep_avg_mprage_N20_mask_3_bilateral_Hippanterior_MNI.nii.gz $anat_dir/"$prefix"HPCanteriorGruber"$suffix"
3dcopy $input_dir/Gruber_HPC/0426rrewrep_avg_mprage_N20_mask_5_bilateral_Hippposterior_MNI.nii.gz $anat_dir/"$prefix"HPCposteriorGruber"$suffix"

### PRIMARY SENSORY CORTICES ###

# use AFNI to extract V1 and A1 mask from Glasser et al 2016 atlas
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:L_Primary_Visual_Cortex -prefix $anat_dir/"$prefix"V1_L"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:R_Primary_Visual_Cortex -prefix $anat_dir/"$prefix"V1_R"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:L_Second_Visual_Area -prefix $anat_dir/"$prefix"V2_L"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:R_Second_Visual_Area -prefix $anat_dir/"$prefix"V2_R"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:L_Primary_Auditory_Cortex -prefix $anat_dir/"$prefix"A1_L"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:R_Primary_Auditory_Cortex -prefix $anat_dir/"$prefix"A1_R"$suffix"


# combine masks for both hemispheres and delete unilateral masks
3dcalc -a $anat_dir/"$prefix"V1_L"$suffix" -b $anat_dir/"$prefix"V1_R"$suffix" -expr 'a+b' -prefix $anat_dir/"$prefix"V1"$suffix"
3dcalc -a $anat_dir/"$prefix"V2_L"$suffix" -b $anat_dir/"$prefix"V2_R"$suffix" -expr 'a+b' -prefix $anat_dir/"$prefix"V2"$suffix"
3dcalc -a $anat_dir/"$prefix"A1_L"$suffix" -b $anat_dir/"$prefix"A1_R"$suffix" -expr 'a+b' -prefix $anat_dir/"$prefix"A1"$suffix"

rm $anat_dir/*_L_* $anat_dir/*_R_*

### OVERLAP FUNCTIONAL REWARD AND ANATOMICAL MASKS ###

# term-based meta-analysis conducted on NeuroSynth.org on 30-01-2020
# output: reward_association-test_z_FDR_0.01 & reward_uniformity-test_z_FDR_0.01

# binarise NeuroSynth masks
#3dmask_tool -input reward_association-test_z_FDR_0.01.nii.gz -prefix reward_association-test.nii.gz 
#3dmask_tool -input reward_uniformity-test_z_FDR_0.01.nii.gz -prefix reward_uniformity-test.nii.gz

# define names
gruber_orig=reward_pFgA_z_FDR_0.05_revinf_329studies.nii.nii.gz
gruber=neurosynth_label-reward_mask.nii.gz
gruber_t1=neurosynth_res-anat_label-reward_mask.nii.gz

3dmask_tool -input $input_dir/$gruber_orig -prefix $ROI_dir/$gruber
# resample to space of MNI (anatomical template)
3dresample -master $template -input $ROI_dir/$gruber -prefix $ROI_dir/$gruber_t1

rm $ROI_dir/$gruber

# arrays with ROI names
anat=(NAcc Caudate VTASN HPC HPCanteriorGruber HPCposteriorGruber V1 V2 A1)
#num_anat=${#anat[@]}
reward_sen="NAcc Caudate VTASN HPC"

# define prefix
prefix_t1=MNI_res-anat_label-
prefix_epi=MNI_res-epi_label-

suffix_reward=_desc-rewardsensitiv_mask.nii.gz
suffix_anat=_desc-anatomical_mask.nii.gz


# for each mask in the anat array
#for (( m=0; m<${num_anat}; m++)); do
for roi in "${anat[@]}"; do


    # for some of the ROIs, determine reward-sensitive areas and resample those and anatomical ROI to t1 and epi grid
	if echo $reward_sen | grep -w $roi > /dev/null; then

	    # resample anatomical ROI to anatomical MNI grid
	    3dresample -master $template -input $anat_dir/"$prefix""$roi""$suffix" -prefix $ROI_dir/"$prefix_t1""$roi""$suffix_anat"

	    # resample anatomical ROI to EPI MNI grid
	    3dresample -master "$ROI_dir"/$epi_ave_mtw -input $anat_dir/"$prefix""$roi""$suffix" -prefix  $ROI_dir/"$prefix_epi""$roi""$suffix_anat"

	    # combine t1-resampled anatomical mask and t1-resampled reward mask (inclusively)
	    3dcalc -a $ROI_dir/"$prefix_t1""$roi""$suffix_anat" -b $ROI_dir/$gruber_t1     \
                -expr 'a*b' -prefix $ROI_dir/"$prefix_t1""$roi""$suffix_reward"

	    # resample inclusive mask to epi grid
	    3dresample -master "$ROI_dir"/$epi_ave_mtw -input $ROI_dir/"$prefix_t1""$roi""$suffix_reward" -prefix $ROI_dir/"$prefix_epi""$roi""$suffix_reward"

    # otherwise simply resample anatomical mask to t1 and epi grid
    else
    
	    # resample anatomical ROI to anatomical MNI grid
	    3dresample -master $template -input $anat_dir/"$prefix""$roi""$suffix" -prefix $ROI_dir/"$prefix_t1""$roi""$suffix"

	    # resample anatomical ROI to EPI MNI grid
	    3dresample -master "$ROI_dir"/$epi_ave_mtw -input $anat_dir/"$prefix""$roi""$suffix" -prefix  $ROI_dir/"$prefix_epi""$roi""$suffix"

    fi
	
done

# create anterior and posterior HPC mask based on recommendations of Poppenk et al (2013)

# use MNI y=-21 to determine uncal apex and hence landmask to divide anterior and posterior HPC
# this translates into -A -117 (adding -117 planes of zero at anterior edge) to create posterior HPC
# and into -P -113 [A-P extent 229 voxels --> 229-117=113] to create anterior HPC

# use zeropad to add zeroes to anterior / posterior part of HPC respectively and to hence create pHPC & aHPC
3dZeropad -A -117 -prefix $ROI_dir/"$prefix_t1"pHPC"$suffix" $ROI_dir/"$prefix_t1"HPC"$suffix_anat"
3dZeropad -P -113  -prefix $ROI_dir/"$prefix_t1"aHPC"$suffix" $ROI_dir/"$prefix_t1"HPC"$suffix_anat"

# resample aHPC & pHPC ROI to EPI MNI grid
3dresample -master "$ROI_dir"/$epi_ave_mtw -input $ROI_dir/"$prefix_t1"aHPC"$suffix" -prefix  $ROI_dir/"$prefix_epi"aHPC"$suffix"
3dresample -master "$ROI_dir"/$epi_ave_mtw -input $ROI_dir/"$prefix_t1"pHPC"$suffix" -prefix  $ROI_dir/"$prefix_epi"pHPC"$suffix"



# as last step combine all ROIs
#3dcalc -a $ROI_dir/HPC_anterior_resampled_EPI.nii.gz -b $ROI_dir/HPC_posterior_resampled_EPI.nii.gz -c $ROI_dir/HPC_resampled_EPI.nii.gz -d $ROI_dir/Caudate_resampled_EPI.nii.gz -e $ROI_dir/NAcc_resampled_EPI.nii.gz -f $ROI_dir/SNVTA_resampled_EPI.nii.gz -g $ROI_dir/V1_resampled_EPI.nii.gz -h $ROI_dir/V2_resampled_EPI.nii.gz -i $ROI_dir/A1_resampled_EPI.nii.gz -j $ROI_dir/GM_mask_MNI152_T1_2009c_dilated_resampled.nii.gz  -expr 'a+b+c+d+e+f+g+h+i+j' -prefix $ROI_dir/wholeBrain.nii.gz
#3dmask_tool -input $ROI_dir/wholeBrain.nii.gz -prefix $ROI_dir/wholeBrain_mask.nii.gz
#rm $ROI_dir/wholeBrain.nii.gz


