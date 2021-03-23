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
ana_root=$deriv_dir/analysis/"$task"
ana_path="$ana_root"/3dISC
CS_path=$ana_root/3dClustSim

script_path=$path/scripts/fmri/bold/analysis/3dISC
script_out=$script_path/RSA

# define data table set
covarset=(memo unique)
covarset=(unique)

###### define masks to loop through ######

# define ROI masks
mask_dir=$deriv_dir/ROI_masks/output

masks=(GM)

mask_files=(MNI_res-epi_label-GM_mask.nii.gz)

num_masks=${#masks[@]}

######### copy ISC MAPS ############

# change directory to ISC folder
cd $ISC_path

# get all files in an array
input=ISC_sub*sub*task-"$task"_z.nii.gz

array=(`find $ISC_path -name "${input}"`)
length_array=${#array[@]}
echo $length_array

# change directory to analyses folder
cd $ana_path

# copy files
for (( f=0; f<${length_array}; f++));
    do

	# define file_id
    file_id=${array[$f]}
	echo $file_id
	# copy file
	#cp $file_id .

done


# note: $ana_path contains all ISC maps and data tables for ISC RSA

# copy 3dClustSim output to directory
cp $CS_path/ClustSim_magictrickwatching* $ana_path/


# loop over masks to create a 3dISC for each mask
for (( m=0; m<${num_masks}; m++)); do

	# define variables
	mask=${masks[$m]}
	mask_f=${mask_files[$m]}


	# run the model for magictrickwatching (no covariates)
	#source "$script_out"/3dISC_"$task"_"$mask" > "$script_out"/output_3dISC_"$task"_"$mask".txt

    #Copy the string 'x' (file) into the dataset(s) giving it the name n
    #3drefit -atrstring AFNI_CLUSTSIM_NN1_1sided file:ClustSim_magictrickwatching.NN1_1sided.niml \
	#-atrstring AFNI_CLUSTSIM_NN2_1sided file:ClustSim_magictrickwatching.NN2_1sided.niml \
	#-atrstring AFNI_CLUSTSIM_NN3_1sided file:ClustSim_magictrickwatching.NN3_1sided.niml \
	#-atrstring AFNI_CLUSTSIM_NN1_2sided file:ClustSim_magictrickwatching.NN1_2sided.niml \
	#-atrstring AFNI_CLUSTSIM_NN2_2sided file:ClustSim_magictrickwatching.NN2_2sided.niml \
	#-atrstring AFNI_CLUSTSIM_NN3_2sided file:ClustSim_magictrickwatching.NN3_2sided.niml \
	#-atrstring AFNI_CLUSTSIM_NN1_bisided file:ClustSim_magictrickwatching.NN1_bisided.niml \
	#-atrstring AFNI_CLUSTSIM_NN2_bisided file:ClustSim_magictrickwatching.NN2_bisided.niml \
	#-atrstring AFNI_CLUSTSIM_NN3_bisided file:ClustSim_magictrickwatching.NN3_bisided.niml \
	#ISC_"$task"_"$mask"+tlrc


    for approach in "${covarset[@]}"; do

        ###### define covariates and interaction to loop through ######
        if [[ "$approach" == *"memo"* ]]; then

            # use covariates and interactions from "default" / memo approach
            covariates=(corr_curiosity corr_confidence corr_cuedRecallStrict curBeta_cuedRecallStrict curBetaOnly_cuedRecallStrict curBen_cuedRecallStrict curCor_cuedRecallStrict corr_cuedRecallLenient curBeta_cuedRecallLenient curBetaOnly_cuedRecallLenient curBen_cuedRecallLenient curCor_cuedRecallLenient corr_allConf curBeta_allConf curBetaOnly_allConf curBen_allConf curCor_allConf corr_highConf curBeta_highConf curBetaOnly_highConf curBen_highConf curCor_highConf corr_aboveAvgConf curBeta_aboveAvgConf curBetaOnly_aboveAvgConf curBen_aboveAvgConf curCor_aboveAvgConf corr_rememberedStrictAboveAvg curBeta_rememberedStrictAboveAvg curBetaOnly_rememberedStrictAboveAvg curBen_rememberedStrictAboveAvg curCor_rememberedStrictAboveAvg corr_rememberedLenientAboveAvg curBeta_rememberedLenientAboveAvg curBetaOnly_rememberedLenientAboveAvg curBen_rememberedLenientAboveAvg curCor_rememberedLenientAboveAvg corr_rememberedStrictHigh curBeta_rememberedStrictHigh curBetaOnly_rememberedStrictHigh curBen_rememberedStrictHigh curCor_rememberedStrictHigh corr_rememberedLenientHigh curBeta_rememberedLenientHigh curBetaOnly_rememberedLenientHigh curBen_rememberedLenientHigh curCor_rememberedLenientHigh)

            covariates=(corr_curiosity corr_confidence corr_aboveAvgConf curBeta_aboveAvgConf)

        else

            # use covariates from "unique" approach
            covariates=(uniqueConf uniqueCur_cuedRecallStrict uniqueMem_cuedRecallStrict uniqueBeta_cuedRecallStrict uniqueBetaOnly_cuedRecallStrict uniqueBen_cuedRecallStrict uniqueCor_cuedRecallStrict uniqueCur_cuedRecallLenient uniqueMem_cuedRecallLenient uniqueBeta_cuedRecallLenient uniqueBetaOnly_cuedRecallLenient uniqueBen_cuedRecallLenient uniqueCor_cuedRecallLenient uniqueCur_allConf uniqueMem_allConf uniqueBeta_allConf uniqueBetaOnly_allConf uniqueBen_allConf uniqueCor_allConf uniqueCur_highConf uniqueMem_highConf uniqueBeta_highConf uniqueBetaOnly_highConf uniqueBen_highConf uniqueCor_highConf uniqueCur_aboveAvgConf uniqueMem_aboveAvgConf uniqueBeta_aboveAvgConf uniqueBetaOnly_aboveAvgConf uniqueBen_aboveAvgConf uniqueCor_aboveAvgConf uniqueCur_rememberedStrictAboveAvg uniqueMem_rememberedStrictAboveAvg uniqueBeta_rememberedStrictAboveAvg uniqueBetaOnly_rememberedStrictAboveAvg uniqueBen_rememberedStrictAboveAvg uniqueCor_rememberedStrictAboveAvg uniqueCur_rememberedLenientAboveAvg uniqueMem_rememberedLenientAboveAvg uniqueBeta_rememberedLenientAboveAvg uniqueBetaOnly_rememberedLenientAboveAvg uniqueBen_rememberedLenientAboveAvg uniqueCor_rememberedLenientAboveAvg uniqueCur_rememberedStrictHigh uniqueMem_rememberedStrictHigh uniqueBeta_rememberedStrictHigh uniqueBetaOnly_rememberedStrictHigh uniqueBen_rememberedStrictHigh uniqueCor_rememberedStrictHigh uniqueCur_rememberedLenientHigh uniqueMem_rememberedLenientHigh uniqueBeta_rememberedLenientHigh uniqueBetaOnly_rememberedLenientHigh uniqueBen_rememberedLenientHigh uniqueCor_rememberedLenientHigh)

            covariates=(uniqueConf) # uniqueCur_aboveAvgConf uniqueMem_aboveAvgConf uniqueBeta_aboveAvgConf)
            covariates=(uniqueCur_aboveAvgConf) # uniqueCur_aboveAvgConf uniqueMem_aboveAvgConf uniqueBeta_aboveAvgConf)

        fi

        num_covariates=${#covariates[@]}


        # loop over covariates to run 3dISC for each of them
	    for (( v=0; v<${num_covariates}; v++)); do
	    
            # define variables
            cov=${covariates[$v]}
		    int=${interaction[$v]}

	        # run the RSA model for each covariate
	        #source "$script_out"/3dISC_"$task"_"$cov"_"$mask" > "$script_out"/output_3dISC_"$task"_"$cov"_"$mask".txt


            #Copy the string 'x' (file) into the dataset(s) giving it the name n
            3drefit -atrstring 'AFNI_CLUSTSIM_NN1_1sided' file:ClustSim_magictrickwatching.NN1_1sided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN2_1sided' file:ClustSim_magictrickwatching.NN2_1sided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN3_1sided' file:ClustSim_magictrickwatching.NN3_1sided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN1_2sided' file:ClustSim_magictrickwatching.NN1_2sided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN2_2sided' file:ClustSim_magictrickwatching.NN2_2sided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN3_2sided' file:ClustSim_magictrickwatching.NN3_2sided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN1_bisided' file:ClustSim_magictrickwatching.NN1_bisided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN2_bisided' file:ClustSim_magictrickwatching.NN2_bisided.niml \
	            -atrstring 'AFNI_CLUSTSIM_NN3_bisided' file:ClustSim_magictrickwatching.NN3_bisided.niml \
	        ISC_"$task"_"$cov"_"$mask"+tlrc



	    done
    done
done

# remove all ISC maps and clust stim output from directory
#rm $ana_path/$input
#rm $ana_path/ClustSim_magictrickwatching*



