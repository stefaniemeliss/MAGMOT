#!/bin/tcsh -xef

set subj = sub-control003_task-rest
set output_dir = $subj.results

cd /storage/shared/research/cinn/2018/MAGMOT/derivatives/afniproc/sub-control003/$output_dir

rm rm.*

# set list of runs
set runs = (`count -digits 2 1 2`)

foreach run ( $runs )
    3dAutomask -clfrac 0.15 -dilate 1 -prefix rm.mask_r$run pb03.$subj.r$run.volreg+tlrc
end

# create union of inputs, output type is byte
3dmask_tool -inputs rm.mask_r*+tlrc.HEAD -union -prefix full_mask_015d.$subj

# ---- create subject anatomy mask, mask_anat.$subj+tlrc ----
#      (resampled from tlrc anat)
3dresample -master full_mask_015d.$subj+tlrc -input anatQQ.sub-control003+tlrc \
           -prefix rm.resam.anat

# convert to binary anat mask; fill gaps and holes
3dmask_tool -dilate_input 5 -5 -fill_holes -input rm.resam.anat+tlrc      \
            -prefix mask_anat_015d.$subj

# compute tighter EPI mask by intersecting with anat mask
3dmask_tool -input full_mask_015d.$subj+tlrc mask_anat_015d.$subj+tlrc              \
            -inter -prefix mask_epi_anat_015d.$subj


#3dresample -rmode NN -master pb03.sub-control003_task-rest.r01.volreg+tlrc. -input anatQQ_mask.nii.gz -prefix anatQQ_res-epi_mask.nii.gz

#3dcalc -a anatQQ.sub-control003+tlrc. -expr 'astep(a,4)' -prefix anatQQ_mask.nii.gz


#3dmask_tool -union -prefix anat_func_mask033.nii.gz -inputs anatQQ_res-epi_mask.nii.gz full_mask_033.sub-control003_task-rest+tlrc.

