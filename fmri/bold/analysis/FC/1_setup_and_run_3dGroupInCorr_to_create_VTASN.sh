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

############### specify masks and p value ###############

# specify mask (created in previous script)
epi_mask=sample_label-dilatedGM_mask.nii.gz
gm_mask=sample_label-gm_mask.nii.gz

# define number of voxel in mask
nvox=$(3dBrickStat -count -non-zero $ROI_dir/$gm_mask)
ntests=1
ntot=$(echo "$nvox * $ntests" | bc -l)

# bonferoni correction
p=0.05
p_cor=$(echo "$p / $ntot " | bc -l)
echo $p_cor


############### create loop for FC analysis for different smoothing kernels ###############

# define flow variables
run_setup=false
run_corr=false
run_data=false
clustsim=false
# define different smoothing kernels
preproc=(s0 s4 s6 s8)
#preproc=(s0)

# define prefix
task="rest"

# define task runs
run1=$task"_run-1"
run2=$task"_run-2"

# define groups
group1=control
group2=experimental


for smooth in "${preproc[@]}"; do

	############### set up the file used to compute correlations (pearson) ###############

    if ($run_setup == true) then

		cd "$deriv_dir" # change directory to where pre-processed files are

		# --------- PRE - WHOLE SAMPLE ---------

		out_GroupInCorr_pre=FC_task-"$run1"_desc-"$smooth"preproc

		# 3dSetupGroupInCorr: Extract data for 3dGroupInCorr 
		3dSetupGroupInCorr -mask $ROI_dir/$epi_mask -prefix "$FC_dir"/$out_GroupInCorr_pre -byte 	\
			-labels "$code_dir"/labels_r1_"$smooth".txt ./sub-*/func/*_task-"$run1"_desc-"$smooth"preproc_bold.nii.gz

		# --------- PRE - PER GROUP ---------

		out_GroupInCorr_cont_pre=FC_group-"$group1"_task-"$run1"_desc-"$smooth"preproc

		# 3dSetupGroupInCorr: Extract data for 3dGroupInCorr 
		3dSetupGroupInCorr -mask $ROI_dir/$epi_mask -prefix "$FC_dir"/$out_GroupInCorr_cont_pre -byte 	\
			-labels "$code_dir"/labels_r1_cont.txt ./sub-"$group1"*/func/*_task-"$run1"_desc-"$smooth"preproc_bold.nii.gz

		out_GroupInCorr_exp_pre=FC_group-"$group2"_task-"$run1"_desc-"$smooth"preproc

		# 3dSetupGroupInCorr: Extract data for 3dGroupInCorr 
		3dSetupGroupInCorr -mask $ROI_dir/$epi_mask -prefix "$FC_dir"/$out_GroupInCorr_exp_pre -byte 	\
			-labels "$code_dir"/labels_r1_exp.txt ./sub-"$group2"*/func/*_task-"$run1"_desc-"$smooth"preproc_bold.nii.gz


		# --------- POST - WHOLE SAMPLE ---------

		out_GroupInCorr_post=FC_task-"$run2"_desc-"$smooth"preproc

		# 3dSetupGroupInCorr: Extract data for 3dGroupInCorr 
		3dSetupGroupInCorr -mask $ROI_dir/$epi_mask -prefix "$FC_dir"/$out_GroupInCorr_post -byte 		\
			-labels "$code_dir"/labels_r2_"$smooth".txt ./sub-*/func/*_task-"$run2"_desc-"$smooth"preproc_bold.nii.gz

		# --------- POST - PER GROUP ---------

		out_GroupInCorr_cont_post=FC_group-"$group1"_task-"$run2"_desc-"$smooth"preproc

		# 3dSetupGroupInCorr: Extract data for 3dGroupInCorr 
		3dSetupGroupInCorr -mask $ROI_dir/$epi_mask -prefix "$FC_dir"/$out_GroupInCorr_cont_post -byte 	\
			-labels "$code_dir"/labels_r2_cont.txt ./sub-"$group1"*/func/*_task-"$run2"_desc-"$smooth"preproc_bold.nii.gz


		out_GroupInCorr_exp_post=FC_group-"$group2"_task-"$run2"_desc-"$smooth"preproc

		# 3dSetupGroupInCorr: Extract data for 3dGroupInCorr 
		3dSetupGroupInCorr -mask $ROI_dir/$epi_mask -prefix "$FC_dir"/$out_GroupInCorr_exp_post -byte 	\
			-labels "$code_dir"/labels_r2_exp.txt ./sub-"$group2"*/func/*_task-"$run2"_desc-"$smooth"preproc_bold.nii.gz

	fi

	############### compute functional connectivity using HPC as seed ###############
	
	cd "$FC_dir"

    if ($run_corr == true) then

		# --------- PAIRED T-TEST - WHOLE SAMPLE ---------		

		3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr_post.grpincorr.niml -labelA "post_all"		\
						-setB "$FC_dir"/$out_GroupInCorr_pre.grpincorr.niml -labelB "pre_all"	\
						-sendall -verb -paired -clust ClustSim_"$task"_"$smooth"				\
						-batch MASKAVE "$code_dir"/seeds_"$smooth"_paired

		# --------- PAIRED T-TEST - CONTROL ---------		

		3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr_cont_post.grpincorr.niml -labelA "post_cont"		\
						-setB "$FC_dir"/$out_GroupInCorr_cont_pre.grpincorr.niml -labelB "pre_cont"	\
						-sendall -verb -paired															\
						-batch MASKAVE "$code_dir"/seeds_"$smooth"_cont

		# --------- PAIRED T-TEST - EXPERIMENTAL ---------		

		3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr_exp_post.grpincorr.niml -labelA "post_exp"		\
						-setB "$FC_dir"/$out_GroupInCorr_exp_pre.grpincorr.niml -labelB "pre_exp"	\
						-sendall -verb -paired														\
						-batch MASKAVE "$code_dir"/seeds_"$smooth"_exp

		# --------- TWO-SAMPLE T-TEST - PRE ---------

		3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr_cont_pre.grpincorr.niml -labelA "pre_cont"	\
						-setB "$FC_dir"/$out_GroupInCorr_exp_pre.grpincorr.niml -labelB "pre_exp"	\
						-sendall -verb -unpooled													\
						-batch MASKAVE "$code_dir"/seeds_"$smooth"_run1

		# --------- TWO-SAMPLE T-TEST - POST ---------

		3dGroupInCorr -setA "$FC_dir"/$out_GroupInCorr_cont_pre.grpincorr.niml -labelA "post_cont"	\
						-setB "$FC_dir"/$out_GroupInCorr_exp_pre.grpincorr.niml -labelB "post_exp"	\
						-sendall -verb -unpooled													\
						-batch MASKAVE "$code_dir"/seeds_"$smooth"_run2


	fi

	############### extract data for each subject ###############

    if ($run_data == true) then

		mkdir "$smooth"_aHPC

		data=RSFC_paired_"$smooth"_aHPC

		end=105

		for i in $(seq 6 $end); do

			# determine name of sub-brick label
			subbrick_label=$(3dinfo -label "$data"+tlrc.[$i])
			echo $subbrick_label

			# extract sub-brick and save as new file
			3dbucket -prefix "$smooth"_aHPC/"$subbrick_label".nii.gz "$data"+tlrc.[$i]

		done

	############### compute change for each subject in HPC seed connectivity ###############


		cd "$smooth"_aHPC

		# create array with one set of files
		post_files=($(ls A_*))
		pre_files=($(ls B_*))

		num_files="${#post_files[@]}"
		echo "${#pre_files[@]}"
		echo "${#post_files[@]}"

		for ((f=0; f<$num_files; f++)); do

			# determine pre- and post file
			pre_file="${pre_files[f]}"
			post_file="${post_files[f]}"

			#echo $pre_file
			#echo $post_file

			# create prefix
			searchstr="A"
			replacestr="diff"
			diff_prefix="${post_file/$searchstr/$replacestr}"

			# create diff file
			if [ ! -f "$diff_prefix" ]; then

				echo $diff_prefix

				# calculate difference
				3dcalc -a $post_file -b $pre_file -expr 'a-b' -prefix $diff_prefix

			fi
		

		done

	fi

	cd $FC_dir


	############### Two Sample t-test for each task seperately and the combined maps ###############


	# run 3dttest++ with Clustsim for thresholding
	if ($clustsim); then 
		3dttest++ -setA "$smooth"_aHPC/diff_c*  -setB "$smooth"_aHPC/diff_e* -labelA $group1 -labelB $group2 -unpooled -mask $ROI_dir/$gm_mask -prefix 2sample_diff_aHPC_"$smooth" -ClustSim
	#else
		#3dttest++ -setA "$smooth"_aHPC/diff_c*  -setB "$smooth"_aHPC/diff_e* -labelA $group1 -labelB $group2 -unpooled -mask $ROI_dir/$gm_mask -prefix 2sample_diff_aHPC_"$smooth"
	fi	


	############### create VTASN mask ###############

    # inclusively mask the cluster map and the anatomical VTA ROI
    VTASN=MNI_res-epi_label-VTASN_desc-anatomical_mask.nii.gz

	# correct for multiple comparisons within midbrain ROI

	# define number of voxel in mask
	vox_vta=$(3dBrickStat -count -non-zero $ROI_dir/$VTASN)
	ntests=1
	ntot=$(echo "$vox_vta * $ntests" | bc -l)

	# bonferoni correction
	p=0.05
	p_vta=$(echo "$p / $ntot " | bc -l)
	echo $p_vta

	# copy clustsim output
	cp $CS_dir/ClustSim_"$task"_"$smooth"*  .

	# define HPC ROIS
	HPC_ROIs=(RSFC_paired_"$smooth"_aHPC RSFC_paired_"$smooth"_rewardHPC)

	num_ROIs=${#HPC_ROIs[@]}

	# loop over ROIs
	for (( r=0; r<${num_ROIs}; r++)); do

		# copy information from 3dClustSim into files
		3drefit -atrstring 'AFNI_CLUSTSIM_NN1_1sided' file:ClustSim_rest_"$smooth".NN1_1sided.niml  \
		-atrstring AFNI_CLUSTSIM_NN2_1sided file:ClustSim_rest_"$smooth".NN2_1sided.niml      \
		-atrstring AFNI_CLUSTSIM_NN3_1sided file:ClustSim_rest_"$smooth".NN3_1sided.niml      \
		-atrstring AFNI_CLUSTSIM_NN1_2sided file:ClustSim_rest_"$smooth".NN1_2sided.niml      \
		-atrstring AFNI_CLUSTSIM_NN2_2sided file:ClustSim_rest_"$smooth".NN2_2sided.niml      \
		-atrstring AFNI_CLUSTSIM_NN3_2sided file:ClustSim_rest_"$smooth".NN3_2sided.niml      \
		-atrstring AFNI_CLUSTSIM_NN1_bisided file:ClustSim_rest_"$smooth".NN1_bisided.niml    \
		-atrstring AFNI_CLUSTSIM_NN2_bisided file:ClustSim_rest_"$smooth".NN2_bisided.niml    \
		-atrstring AFNI_CLUSTSIM_NN3_bisided file:ClustSim_rest_"$smooth".NN3_bisided.niml    \
		"${HPC_ROIs[$r]}"+tlrc


		# output 3dinfor -label "${HPC_ROIs[$r]}"
		# 0 post_all-pre_all_mean
		#1 post_all-pre_all_Zscr
		#2 post_all_means
		#3 post_all_Zscr
		#4 pre_all_mean
		#5 pre_all_Zscr
	
		# Alpha --> Cluster threshold for $p_vta
		# s0: 0.10 -> 16; 0.05 -> 20; 0.01 -> 29
		# s4: 0.10 -> 16; 0.05 -> 20; 0.01 -> 28
		# s6: 0.10 -> 16; 0.05 -> 20; 0.01 -> 28
		# s8: 0.10 -> 16; 0.05 -> 20; 0.01 -> 28

		
		# create a map with voxels showing significant FC with bilateral HPC at pre-learning rest within the VTA mask (surviving bonferroni correction)
    	#3dClusterize -inset "${HPC_ROIs[$r]}"+tlrc -ithr 5 -idat 4 -NN 1 -clust_nvox 20 -bisided p=$p_vta -pref_map VTA_${HPC_ROIs[$r]}.nii.gz -mask $ROI_dir/$VTASN -binary
    	3dClusterize -inset "${HPC_ROIs[$r]}"+tlrc -ithr 5 -idat 4 -NN 1 -bisided p=$p_vta -pref_map VTA_${HPC_ROIs[$r]}.nii.gz -mask $ROI_dir/$VTASN -binary


	done

	#rm $FC_dir/ClustSim_"$task"_"$smooth"*

done



#####  for BNA POSTER

# extract map
3dClusterize -inset 2sample_diff_aHPC_s4+tlrc -ithr 3 -idat 2 -NN 1 -bisided p=0.001 -clust_nvox 7 -pref_dat change_hpc_fc_control.nii.gz -mask $ROI_dir/$gm_mask -binary

# create image (FIGURE 5)
@chauffeur_afni -ulay $template -olay RSFC_paired_s4_aHPC+tlrc -thr_olay_p2stat $p_cor -thr_olay_pside 2sided -set_subbricks 0 4 5 -cbar Spectrum:red_to_blue -pbar_posonly -prefix aHPC_s4_pre -pbar_saveim aHPC_s4_pre.png -pbar_dim 64x1351H -func_range 0.4 -opacity 6 -label_mode 1 -label_setback 0.5 -label_color black -zerocolor white -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY

@chauffeur_afni -ulay $template -olay RSFC_paired_s4_aHPC+tlrc -thr_olay_p2stat $p_cor -thr_olay_pside 2sided -set_subbricks 0 2 3 -cbar Spectrum:red_to_blue -pbar_posonly -prefix aHPC_s4_post -pbar_saveim aHPC_s4_post.png -pbar_dim 64x1351H -func_range 0.4 -opacity 6 -label_mode 1 -label_setback 0.5 -label_color black -zerocolor white -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY

@chauffeur_afni -ulay $template -olay RSFC_run-2_s4_aHPC+tlrc -thr_olay_p2stat $p_cor -thr_olay_pside 2sided -set_subbricks 0 2 3 -cbar Spectrum:red_to_blue -pbar_posonly -prefix aHPC_s4_post_c -pbar_saveim aHPC_s4_post_c.png -pbar_dim 64x1351H -func_range 0.4 -opacity 6 -label_mode 1 -label_setback 0.5 -label_color black -zerocolor white -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY


@chauffeur_afni -ulay $template -olay RSFC_run-2_s4_aHPC+tlrc -thr_olay_p2stat $p_cor -thr_olay_pside 2sided -set_subbricks 0 4 5 -cbar Spectrum:red_to_blue -pbar_posonly -prefix aHPC_s4_post_e -pbar_saveim aHPC_s4_post_e.png -pbar_dim 64x1351H -func_range 0.4 -opacity 6 -label_mode 1 -label_setback 0.5 -label_color black -zerocolor white -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY

# show effects of smoothing (included in BNA poster presentation slides)
@chauffeur_afni -ulay $template -olay RSFC_paired_s0_aHPC+tlrc -thr_olay_p2stat $p_vta -thr_olay_pside 2sided -set_subbricks 0 4 5 -cbar Spectrum:red_to_blue -pbar_posonly -prefix sFC_s0_pre -pbar_saveim sFC_s0_pre.png -pbar_dim 64x1351H -func_range 0.3 -opacity 6 -label_mode 1 -label_setback 0.5 -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY

@chauffeur_afni -ulay $template -olay RSFC_paired_s4_aHPC+tlrc -thr_olay_p2stat $p_vta -thr_olay_pside 2sided -set_subbricks 0 4 5 -cbar Spectrum:red_to_blue -pbar_posonly -prefix sFC_s4_pre -pbar_saveim sFC_s4_pre.png -pbar_dim 64x1351H -func_range 0.3 -opacity 6 -label_mode 1 -label_setback 0.5 -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY

@chauffeur_afni -ulay $template -olay RSFC_paired_s6_aHPC+tlrc -thr_olay_p2stat $p_vta -thr_olay_pside 2sided -set_subbricks 0 4 5 -cbar Spectrum:red_to_blue -pbar_posonly -prefix sFC_s6_pre -pbar_saveim sFC_s6_pre.png -pbar_dim 64x1351H -func_range 0.3 -opacity 6 -label_mode 1 -label_setback 0.5 -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY

@chauffeur_afni -ulay $template -olay RSFC_paired_s8_aHPC+tlrc -thr_olay_p2stat $p_vta -thr_olay_pside 2sided -set_subbricks 0 4 5 -cbar Spectrum:red_to_blue -pbar_posonly -prefix sFC_s8_pre -pbar_saveim sFC_s8_pre.png -pbar_dim 64x1351H -func_range 0.3 -opacity 6 -label_mode 1 -label_setback 0.5 -montx 7 -monty 2 -box_focus_slices AMASK_FOCUS_OLAY


