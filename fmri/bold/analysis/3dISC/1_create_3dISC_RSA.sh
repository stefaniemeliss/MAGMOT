#!/bin/bash

module load afni19.3.03

source ~/.bashrc

# this script creates the 3dISC command the ISC-RSA analyses

# define task
task=magictrickwatching

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"
deriv_dir="$path"/derivatives
ISC_path="$deriv_dir"/ISC
ana_path="$deriv_dir"/analysis/$task/3dISC
mkdir $ana_path

script_path=$path/scripts/fmri/bold/analysis/3dISC
script_out=$script_path/RSA
mkdir -p $script_out
cd $script_out


# define data table set
covarset=(memo unique)


###### define masks to loop through ######

# define ROI masks
mask_dir=$deriv_dir/ROI_masks/output

#-masks=(pHPC aHPC HPC VTA Caudate NAcc GM)
#masks=(V1 V2 A1 Caudate_func)
masks=(GM)

#mask_files=(HPC_posterior_resampled_EPI.nii.gz HPC_anterior_resampled_EPI.nii.gz HPC_resampled_EPI.nii.gz SNVTA_resampled_EPI.nii.gz Caudate_resampled_EPI.nii.gz NAcc_resampled_EPI.nii.gz wholeBrain_mask.nii.gz)
# define ROI masks: visual, auditory, Nucleus Caudate
#mask_files=(V1_resampled_EPI.nii.gz V2_resampled_EPI.nii.gz A1_resampled_EPI.nii.gz Caudate_func_resampled_EPI.nii.gz)

mask_files=(MNI_res-epi_label-GM_mask.nii.gz)

# redo with new ROI masks
#masks=(GruberVTA GruberNAcc GruberHPC)
#mask_files=(VTASN_Gruber_resampled.nii.gz NAcc_Gruber_resampled.nii.gz HPC_Gruber_resampled.nii.gz)

num_masks=${#masks[@]}




# loop over masks to create a 3dISC for each mask
for (( m=0; m<${num_masks}; m++)); do

	# define variables
	mask=${masks[$m]}
	mask_f=${mask_files[$m]}

	# specify 3dISC command
	printf "3dISC -prefix ISC_"$task"_"$mask" -jobs 4 \\" > ./3dISC_"$task"_"$mask" #name and number of jobs
	printf "\n\t -model 'grp+(1|Subj1)+(1|Subj2)' \\" >> ./3dISC_"$task"_"$mask" # model
	printf "\n\t -qVars   'grp' \\" >> ./3dISC_"$task"_"$mask" # qVars
	printf "\n\t -gltCode ave '1 0' \\" >> ./3dISC_"$task"_"$mask" # contrasts
	printf "\n\t -gltCode G11vG22 '0 1' \\" >> ./3dISC_"$task"_"$mask" # contrasts
	printf "\n\t -gltCode G11 '1 0.5' \\" >> ./3dISC_"$task"_"$mask" # contrasts
	printf "\n\t -gltCode G22 '1 -0.5' \\" >> ./3dISC_"$task"_"$mask" # contrasts
	printf "\n\t -mask $mask_dir$mask_f \\" >> ./3dISC_"$task"_"$mask" # mask
	printf "\n\t -dataTable @dataTable_magictrickwatching_memo.txt \\" >> ./3dISC_"$task"_"$mask" # datatable

    for approach in "${covarset[@]}"; do

        ###### define covariates and interaction to loop through ######
        if [[ "$approach" == *"memo"* ]]; then

            # use covariates and interactions from "default" / memo approach
            covariates=(corr_curiosity corr_confidence corr_cuedRecallStrict curBeta_cuedRecallStrict curBetaOnly_cuedRecallStrict curBen_cuedRecallStrict curCor_cuedRecallStrict corr_cuedRecallLenient curBeta_cuedRecallLenient curBetaOnly_cuedRecallLenient curBen_cuedRecallLenient curCor_cuedRecallLenient corr_allConf curBeta_allConf curBetaOnly_allConf curBen_allConf curCor_allConf corr_highConf curBeta_highConf curBetaOnly_highConf curBen_highConf curCor_highConf corr_aboveAvgConf curBeta_aboveAvgConf curBetaOnly_aboveAvgConf curBen_aboveAvgConf curCor_aboveAvgConf corr_rememberedStrictAboveAvg curBeta_rememberedStrictAboveAvg curBetaOnly_rememberedStrictAboveAvg curBen_rememberedStrictAboveAvg curCor_rememberedStrictAboveAvg corr_rememberedLenientAboveAvg curBeta_rememberedLenientAboveAvg curBetaOnly_rememberedLenientAboveAvg curBen_rememberedLenientAboveAvg curCor_rememberedLenientAboveAvg corr_rememberedStrictHigh curBeta_rememberedStrictHigh curBetaOnly_rememberedStrictHigh curBen_rememberedStrictHigh curCor_rememberedStrictHigh corr_rememberedLenientHigh curBeta_rememberedLenientHigh curBetaOnly_rememberedLenientHigh curBen_rememberedLenientHigh curCor_rememberedLenientHigh)

            interaction=(grCorr_curiosity grCorr_confidence grCorr_cuedRecallStrict grCurBeta_cuedRecallStrict grCurBetaOnly_cuedRecallStrict grCurBen_cuedRecallStrict grCurCor_cuedRecallStrict grCorr_cuedRecallLenient grCurBeta_cuedRecallLenient grCurBetaOnly_cuedRecallLenient grCurBen_cuedRecallLenient grCurCor_cuedRecallLenient grCorr_allConf grCurBeta_allConf grCurBetaOnly_allConf grCurBen_allConf grCurCor_allConf grCorr_highConf grCurBeta_highConf grCurBetaOnly_highConf grCurBen_highConf grCurCor_highConf grCorr_aboveAvgConf grCurBeta_aboveAvgConf grCurBetaOnly_aboveAvgConf grCurBen_aboveAvgConf grCurCor_aboveAvgConf grCorr_rememberedStrictAboveAvg grCurBeta_rememberedStrictAboveAvg grCurBetaOnly_rememberedStrictAboveAvg grCurBen_rememberedStrictAboveAvg grCurCor_rememberedStrictAboveAvg grCorr_rememberedLenientAboveAvg grCurBeta_rememberedLenientAboveAvg grCurBetaOnly_rememberedLenientAboveAvg grCurBen_rememberedLenientAboveAvg grCurCor_rememberedLenientAboveAvg grCorr_rememberedStrictHigh grCurBeta_rememberedStrictHigh grCurBetaOnly_rememberedStrictHigh grCurBen_rememberedStrictHigh grCurCor_rememberedStrictHigh grCorr_rememberedLenientHigh grCurBeta_rememberedLenientHigh grCurBetaOnly_rememberedLenientHigh grCurBen_rememberedLenientHigh grCurCor_rememberedLenientHigh)

            data_table=@dataTable_magictrickwatching_memo.txt

        else

            # use covariates and interactions from "unique" approach
            covariates=(uniqueConf uniqueCur_cuedRecallStrict uniqueMem_cuedRecallStrict uniqueBeta_cuedRecallStrict uniqueBetaOnly_cuedRecallStrict uniqueBen_cuedRecallStrict uniqueCor_cuedRecallStrict uniqueCur_cuedRecallLenient uniqueMem_cuedRecallLenient uniqueBeta_cuedRecallLenient uniqueBetaOnly_cuedRecallLenient uniqueBen_cuedRecallLenient uniqueCor_cuedRecallLenient uniqueCur_allConf uniqueMem_allConf uniqueBeta_allConf uniqueBetaOnly_allConf uniqueBen_allConf uniqueCor_allConf uniqueCur_highConf uniqueMem_highConf uniqueBeta_highConf uniqueBetaOnly_highConf uniqueBen_highConf uniqueCor_highConf uniqueCur_aboveAvgConf uniqueMem_aboveAvgConf uniqueBeta_aboveAvgConf uniqueBetaOnly_aboveAvgConf uniqueBen_aboveAvgConf uniqueCor_aboveAvgConf uniqueCur_rememberedStrictAboveAvg uniqueMem_rememberedStrictAboveAvg uniqueBeta_rememberedStrictAboveAvg uniqueBetaOnly_rememberedStrictAboveAvg uniqueBen_rememberedStrictAboveAvg uniqueCor_rememberedStrictAboveAvg uniqueCur_rememberedLenientAboveAvg uniqueMem_rememberedLenientAboveAvg uniqueBeta_rememberedLenientAboveAvg uniqueBetaOnly_rememberedLenientAboveAvg uniqueBen_rememberedLenientAboveAvg uniqueCor_rememberedLenientAboveAvg uniqueCur_rememberedStrictHigh uniqueMem_rememberedStrictHigh uniqueBeta_rememberedStrictHigh uniqueBetaOnly_rememberedStrictHigh uniqueBen_rememberedStrictHigh uniqueCor_rememberedStrictHigh uniqueCur_rememberedLenientHigh uniqueMem_rememberedLenientHigh uniqueBeta_rememberedLenientHigh uniqueBetaOnly_rememberedLenientHigh uniqueBen_rememberedLenientHigh uniqueCor_rememberedLenientHigh)

            interaction=(grUniqueConf grUniqueCur_cuedRecallStrict grUniqueMem_cuedRecallStrict grUniqueBeta_cuedRecallStrict grUniqueBetaOnly_cuedRecallStrict grUniqueBen_cuedRecallStrict grUniqueCor_cuedRecallStrict grUniqueCur_cuedRecallLenient grUniqueMem_cuedRecallLenient grUniqueBeta_cuedRecallLenient grUniqueBetaOnly_cuedRecallLenient grUniqueBen_cuedRecallLenient grUniqueCor_cuedRecallLenient grUniqueCur_allConf grUniqueMem_allConf grUniqueBeta_allConf grUniqueBetaOnly_allConf grUniqueBen_allConf grUniqueCor_allConf grUniqueCur_highConf grUniqueMem_highConf grUniqueBeta_highConf grUniqueBetaOnly_highConf grUniqueBen_highConf grUniqueCor_highConf grUniqueCur_aboveAvgConf grUniqueMem_aboveAvgConf grUniqueBeta_aboveAvgConf grUniqueBetaOnly_aboveAvgConf grUniqueBen_aboveAvgConf grUniqueCor_aboveAvgConf grUniqueCur_rememberedStrictAboveAvg grUniqueMem_rememberedStrictAboveAvg grUniqueBeta_rememberedStrictAboveAvg grUniqueBetaOnly_rememberedStrictAboveAvg grUniqueBen_rememberedStrictAboveAvg grUniqueCor_rememberedStrictAboveAvg grUniqueCur_rememberedLenientAboveAvg grUniqueMem_rememberedLenientAboveAvg grUniqueBeta_rememberedLenientAboveAvg grUniqueBetaOnly_rememberedLenientAboveAvg grUniqueBen_rememberedLenientAboveAvg grUniqueCor_rememberedLenientAboveAvg grUniqueCur_rememberedStrictHigh grUniqueMem_rememberedStrictHigh grUniqueBeta_rememberedStrictHigh grUniqueBetaOnly_rememberedStrictHigh grUniqueBen_rememberedStrictHigh grUniqueCor_rememberedStrictHigh grUniqueCur_rememberedLenientHigh grUniqueMem_rememberedLenientHigh grUniqueBeta_rememberedLenientHigh grUniqueBetaOnly_rememberedLenientHigh grUniqueBen_rememberedLenientHigh grUniqueCor_rememberedLenientHigh)

            data_table=@dataTable_magictrickwatching_unique.txt

        fi

        num_covariates=${#covariates[@]}
        num_interaction=${#interaction[@]}


        # loop over covariates to create a 3dISC for each of them
	    for (( v=0; v<${num_covariates}; v++)); do
	    
            # define variables
            cov=${covariates[$v]}
		    int=${interaction[$v]}

		    # specify 3dISC command including group
		    printf "3dISC -prefix ISC_"$task"_"$cov"_"$mask" -jobs 4 \\" > ./3dISC_"$task"_"$cov"_"$mask" #name and number of jobs
		    printf "\n\t -model 'grp+$cov+$int+(1|Subj1)+(1|Subj2)' \\" >> ./3dISC_"$task"_"$cov"_"$mask" # model
		    printf "\n\t -qVars   'grp,$cov,$int' \\" >> ./3dISC_"$task"_"$cov"_"$mask" # qVars
		    printf "\n\t -qVarCenters   '$cov,$int' \\" >> ./3dISC_"$task"_"$cov"_"$mask" # qVarCenters
		    printf "\n\t -gltCode ave '1 0 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -gltCode G11vG22 '0 1 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -gltCode G11 '1 0.5 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -gltCode G22 '1 -0.5 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -gltCode $cov '0 0 1 0'  \\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -gltCode "$cov"1v"$cov"2 '0 0 0 1'\\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -gltCode "$cov"1 '0 0 1 0.5'  \\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -gltCode "$cov"2 '0 0 1 -0.5'  \\" >> ./3dISC_"$task"_"$cov"_"$mask" # contrasts
		    printf "\n\t -mask $mask_dir/$mask_f \\" >> ./3dISC_"$task"_"$cov"_"$mask" # mask
		    printf "\n\t -dataTable $data_table \\" >> ./3dISC_"$task"_"$cov"_"$mask" # datatable

	    done
    done
done






