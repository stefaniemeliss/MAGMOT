#!/bin/tcsh -xef


# Set top level directory structure
set topdir = /storage/shared/research/cinn/2018/MAGMOT #study folder
echo $topdir
#set task = rest
set task = magictrickwatching
set derivroot = $topdir/derivatives
set fsroot = $derivroot/FreeSurfers
set outroot = $derivroot/afniproc

# define subject listecho $
set BIDSdir = $topdir/MAGMOT_BIDS

cd $BIDSdir
set subjects	=(`ls -d sub*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-control*`) # this creates an array containing all subjects in the BIDS directory
set subjects	=(`ls -d sub-experimental*`) # this creates an array containing all subjects in the BIDS directory
echo $subjects
echo $#subjects

set subjects	= sub-experimental016
#set subjects	= sub-control001



# for each subject in the subjects array
foreach subject ($subjects)

    set subj = "$subject"_task-"$task"
    set output_dir = $subj.results

    cd $outroot/$subject/$output_dir

    # set list of runs
    #set runs = (`count -digits 2 1 2`) # rest
    set runs = (`count -digits 2 1 3`) # task
    set runs = (`count -digits 2 1 4`) # sub-experimental016

    foreach run ( $runs )
        3dAutomask -clfrac 0.1 -prefix rm.mask_r$run pb03.$subj.r$run.volreg+tlrc
    end

    # create union of inputs, output type is byte
    3dmask_tool -inputs rm.mask_r*+tlrc.HEAD -union -prefix full_mask_010.$subj

    # ---- create subject anatomy mask, mask_anat.$subj+tlrc ----
    #      (resampled from tlrc anat)
    3dresample -master full_mask_010.$subj+tlrc -input anatQQ."$subject"+tlrc \
               -prefix rm.resam.anat

    # convert to binary anat mask; fill gaps and holes
    3dmask_tool -dilate_input 5 -5 -fill_holes -input rm.resam.anat+tlrc      \
                -prefix mask_anat_010.$subj

    # compute tighter EPI mask by intersecting with anat mask
    3dmask_tool -input full_mask_010.$subj+tlrc mask_anat_010.$subj+tlrc              \
                -inter -prefix mask_epi_anat_010.$subj


    # blur
    cp 3dFWHMx.1D 3dFWHMx_orig.1D
    rm 3dFWHMx.1D

    foreach run ( $runs )
        3dBlurToFWHM -FWHM 8 -mask mask_epi_anat_010.$subj+tlrc                   \
                     -input pb03.$subj.r$run.volreg+tlrc \
                     -prefix pb04.$subj.r$run.blur_010 
    end

    # run the regression analysis
    3dDeconvolve -input pb04.$subj.r*.blur_010+tlrc.HEAD                          \
        -censor censor_${subj}_combined_2.1D                                  \
        -ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
        -ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
        -ortvec mot_demean.r01.1D mot_demean_r01                              \
        -ortvec mot_demean.r02.1D mot_demean_r02                              \
        -ortvec mot_deriv.r01.1D mot_deriv_r01                                \
        -ortvec mot_deriv.r02.1D mot_deriv_r02                                \
        -polort 6                                                             \
        -num_stimts 0                                                         \
        -fout -tout -x1D X_010.xmat.1D -xjpeg X_010.jpg                               \
        -x1D_uncensored X_010.nocensor.xmat.1D                                    \
        -fitts fitts_010.$subj                                                    \
        -errts errts_010.${subj}                                                  \
        -x1D_stop                                                             \
        -bucket stats_010.$subj

    # -- use 3dTproject to project out regression matrix --
    #    (make errts like 3dDeconvolve, but more quickly)
    3dTproject -polort 0 -input pb04.$subj.r*.blur_010+tlrc.HEAD                  \
               -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
               -ort X_010.nocensor.xmat.1D -prefix errts_010.${subj}.tproject

    # --------------------------------------------------
    # generate fast ANATICOR result: errts.$subj.fanaticor+tlrc


    # -- use 3dTproject to project out regression matrix --
    #    (make errts like 3dDeconvolve, but more quickly)
    3dTproject -polort 0 -input pb04.$subj.r*.blur_010+tlrc.HEAD                  \
               -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
               -dsort Local_FSWMe_rall+tlrc                                   \
               -ort X_010.nocensor.xmat.1D -prefix errts_010.$subj.fanaticor

end
