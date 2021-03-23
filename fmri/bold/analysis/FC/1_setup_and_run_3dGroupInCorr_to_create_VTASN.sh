#!/bin/bash
source ~/.bashrc

module load afni19.3.03

# this code runs a whole brain seed-based functional connectivity analysis using the HPC as seed to determine areas that exhibit functional connectivity

# the whole brain map is then intersecte with a VTA map to determine which voxel in the VTA show a significant FC with the HPC at pre learning rest


# define path
path="/storage/shared/research/cinn/2018/MAGMOT"

# define task
task=rest

# define directories
deriv_dir=$path/derivatives
ROI_dir="$deriv_dir"/ROI_masks/output
anal_dir=$deriv_dir/analysis/"$task"
CS_dir=$anal_dir/3dClustSim
FC_dir=$anal_dir/FC
mkdir $FC_dir
code_dir=$path/scripts/fmri/bold/analysis/FC

# 3dSetupGroupInCorr: Extract data for 3dGroupInCorr

cd "$deriv_dir" # change directory to where pre-processed files are

############### set up the file used to compute correlations (pearson) ###############

# define EPI extent mask
union_mask=sample_task-"$task"_label-GM_desc-MNIunion_mask.nii.gz
inter_mask=sample_task-"$task"_label-GM_desc-MNIinter_mask.nii.gz

# define prefix
out_GroupInCorr=FC_task-"$task"_run-1_wholebrain_pearson
out_GroupInCorr_s0=FC_nosmooth_task-"$task"_run-1_wholebrain_pearson
out_GroupInCorr_s4=FC_smoothed_task-"$task"_run-1_wholebrain_pearson

# set up FC using *smoothed* data
#3dSetupGroupInCorr -mask "$ROI_dir"/$inter_mask -prefix "$FC_dir"/$out_GroupInCorr_s4 -byte ./sub*/func/*_task-"$task"_run-1_desc-smoothed_bold.nii.gz
3dSetupGroupInCorr -mask "$ROI_dir"/$inter_mask -prefix "$FC_dir"/$out_GroupInCorr_s4 -byte ./sub*/func/*_task-"$task"_run-1_desc-preproc_bold.nii.gz

# set up FC using *unsmoothed* data 
3dSetupGroupInCorr -mask "$ROI_dir"/$inter_mask -prefix "$FC_dir"/$out_GroupInCorr_s0 -byte ./sub*/func/*_task-"$task"_run-1_desc-nosmooth_bold.nii.gz



#3dSetupGroupInCorr -mask "$ROI_dir"/$inter_mask -prefix "$FC_dir"/$out_GroupInCorr -byte ./sub*/func/*_task-"$task"_run-1_desc-preproc_bold.nii.gz
# set up the file used to compute correlations (spearman)
#3dSetupGroupInCorr -mask "$path"/derivatives/ROI_masks/epi_mask_MAGMOT_restingstate.nii.gz -prefix "$path"/derivatives/restingstate/analysis/RSFC_run-1/RSFC_run-1_wholebrain_spearman -byte -prep SPEARMAN ./*/*_task-restingstate_run-1_afniproc.nii.gz

############### compute functional connectivity using HPC as seed ###############

### At this point you could run (in 2 separate terminal windows)
cd "$FC_dir"

# compute whole brain connectivity with HPC seed using smoothed data
3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr_s4.grpincorr.niml -verb -batch MASKAVE "$code_dir"/seeds_smoothed

# compute whole brain connectivity with HPC seed using smoothed data
3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr_s0.grpincorr.niml -verb -batch MASKAVE "$code_dir"/seeds_nosmooth



# compute whole brain correlation with seed region using Pearson
#3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr.grpincorr.niml -verb -batch MASKAVE "$code_dir"/seeds_pearson
#aHPC_f_RSFC_pearson_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_anterior_func_resampled_EPI.nii.gz
#aHPC_a_RSFC_pearson_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_anterior_resampled_EPI.nii.gz
#HPC_f_RSFC_pearson_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_func_resampled_EPI.nii.gz
#HPC_a_RSFC_pearson_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_resampled_EPI.nii.gz
#pHPC_f_RSFC_pearson_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_posterior_func_resampled_EPI.nii.gz
#pHPC_a_RSFC_pearson_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_posterior_resampled_EPI.nii.gz

# compute whole brain correlation with seed region using Spearman
#3dGroupInCorr -setA "$path"/derivatives/restingstate/analysis/RSFC_run-1/RSFC_run-1_wholebrain_spearman.grpincorr.niml -verb -batch MASKAVE "$path"/derivatives/restingstate/scripts/seeds_spearman
#aHPC_f_RSFC_spearman_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_anterior_func_resampled_EPI.nii.gz
#aHPC_a_RSFC_spearman_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_anterior_resampled_EPI.nii.gz
#HPC_f_RSFC_spearman_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_func_resampled_EPI.nii.gz
#HPC_a_RSFC_spearman_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_resampled_EPI.nii.gz
#pHPC_f_RSFC_spearman_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_posterior_func_resampled_EPI.nii.gz
#pHPC_a_RSFC_spearman_run-1 /storage/shared/research/cinn/2018/MAGMOT/derivatives/ROI_masks/HPC_posterior_resampled_EPI.nii.gz

HPC_ROIs=(HPC_f_RSFC_pearson_run-1 aHPC_a_RSFC_pearson_run-1 HPC_f_RSFC_spearman_run-1 aHPC_a_RSFC_spearman_run-1)
HPC_ROIs=(RSFC_run-1_pearson_aHPC RSFC_run-1_pearson_rewardHPC)
HPC_ROIs=(RSFC_smoothed_run-1_pearson_rewardHPC RSFC_nosmooth_run-1_pearson_rewardHPC RSFC_smoothed_run-1_pearson_aHPC RSFC_nosmooth_run-1_pearson_aHPC)
num_ROIs=${#HPC_ROIs[@]}

for (( r=0; r<${num_ROIs}; r++)); do

    # copy ClustSim output
    cp $CS_dir/ClustSim_rest* $FC_dir/

    # copy information from 3dClustSim into files
    3drefit -atrstring 'AFNI_CLUSTSIM_NN1_1sided' file:ClustSim_rest.NN1_1sided.niml  \
    -atrstring AFNI_CLUSTSIM_NN2_1sided file:ClustSim_rest.NN2_1sided.niml      \
    -atrstring AFNI_CLUSTSIM_NN3_1sided file:ClustSim_rest.NN3_1sided.niml      \
    -atrstring AFNI_CLUSTSIM_NN1_2sided file:ClustSim_rest.NN1_2sided.niml      \
    -atrstring AFNI_CLUSTSIM_NN2_2sided file:ClustSim_rest.NN2_2sided.niml      \
    -atrstring AFNI_CLUSTSIM_NN3_2sided file:ClustSim_rest.NN3_2sided.niml      \
    -atrstring AFNI_CLUSTSIM_NN1_bisided file:ClustSim_rest.NN1_bisided.niml    \
    -atrstring AFNI_CLUSTSIM_NN2_bisided file:ClustSim_rest.NN2_bisided.niml    \
    -atrstring AFNI_CLUSTSIM_NN3_bisided file:ClustSim_rest.NN3_bisided.niml    \
    ${HPC_ROIs[$r]}+tlrc

    rm $FC_dir/ClustSim_rest*

    # save outputs in cluster map with areas showing a significant RSFC with HPC ROI
    if [[ "${HPC_ROIs[$r]}" == *"smoothed"* ]]; then
        # Use second-nearest neighbor clustering: voxels cluster together if faces OR edges touch
        # bisided: positive and negative voxels above threshold are clustered separately but both sets of clusters are kept
        # 33 voxels --> clustsim output for p = 0.001 and alpha = 0.05
    
        3dClusterize -inset "${HPC_ROIs[$r]}"+tlrc -ithr 1 -idat 0 -NN 2 -bisided p=0.001 -clust_nvox 33 -pref_map Cluster_${HPC_ROIs[$r]}.nii.gz

    else

        3dClusterize -inset "${HPC_ROIs[$r]}"+tlrc -ithr 1 -idat 0 -NN 2 -bisided p=0.001 -clust_nvox 5 -pref_map Cluster_${HPC_ROIs[$r]}.nii.gz

    fi


    # inclusively mask the cluster map and the anatomical VTA ROI
    VTASN=MNI_res-epi_label-VTASN_desc-anatomical_mask.nii.gz

    3dcalc -a Cluster_${HPC_ROIs[$r]}.nii.gz -b $ROI_dir/$VTASN -expr 'a*b' -prefix VTA_${HPC_ROIs[$r]}.nii.gz

done




