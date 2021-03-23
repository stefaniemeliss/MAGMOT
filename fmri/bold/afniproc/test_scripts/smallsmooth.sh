#!/bin/tcsh -xef


# Set top level directory structure
set topdir = /storage/shared/research/cinn/2018/MAGMOT #study folder
echo $topdir
#set task = rest
#set task = magictrickwatching
set tasks = (rest magictrickwatching)
set tasks = (rest)
#set tasks = (magictrickwatching)
set derivroot = $topdir/derivatives
set fsroot = $derivroot/FreeSurfers
set outroot = $derivroot/afniproc

# define subject listecho $
set BIDSdir = $topdir/MAGMOT_BIDS

cd $BIDSdir
set subjects	=(`ls -d sub*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-control*`) # this creates an array containing all subjects in the BIDS directory
#set subjects	=(`ls -d sub-experimental*`) # this creates an array containing all subjects in the BIDS directory
echo $subjects
echo $#subjects

#set subjects	= (sub-experimental050)

#set subjects	= (sub-control001)


foreach task ($tasks)

    # for each subject in the subjects array
    foreach subject ($subjects)

        set subj = "$subject"_task-"$task"
        set output_dir = $subj.results_old

        cd $outroot/$subject/$output_dir

        # set list of runs
        if ( $task == magictrickwatching ) then
            set runs = (`count -digits 2 1 3`) # task
            if ( $subject == sub-experimental016 ) then
                set runs = (`count -digits 2 1 4`) # sub-experimental016
            endif
        else 
            set runs = (`count -digits 2 1 2`) # rest
        endif
        

    foreach run ( $runs )
#        3dBlurToFWHM -FWHM 8 -mask mask_epi_anat_010.$subj+tlrc                   \
#                     -input pb03.$subj.r$run.volreg+tlrc \
#                     -prefix pb04.$subj.r$run.blur

    3dmerge -1blur_fwhm 2 -doall -prefix pb04.$subj.r$run.blursmall \
            pb03.$subj.r$run.volreg+tlrc

    end

    # run the regression analysis
    3dDeconvolve -input pb04.$subj.r*.blursmall+tlrc.HEAD                          \
            -mask mask_epi_anat_010.$subj+tlrc                              \
            -censor censor_${subj}_combined_2.1D                                  \
            -ortvec ROIPC.FSvent.r01.1D ROIPC.FSvent.r01                          \
            -ortvec ROIPC.FSvent.r02.1D ROIPC.FSvent.r02                          \
            -ortvec mot_demean.r01.1D mot_demean_r01                              \
            -ortvec mot_demean.r02.1D mot_demean_r02                              \
            -ortvec mot_deriv.r01.1D mot_deriv_r01                                \
            -ortvec mot_deriv.r02.1D mot_deriv_r02                                \
            -polort 6                                                             \
            -num_stimts 0                                                         \
            -fout -tout -x1D X_smoothed.xmat.1D -xjpeg X_smoothed.jpg                               \
            -x1D_uncensored X_smoothed.nocensor.xmat.1D                                    \
            -fitts fitts_smoothed.$subj                                                    \
            -errts errts_smoothed.${subj}                                                  \
            -x1D_stop                                                             \
            -bucket stats_smoothed.$subj

        # -- use 3dTproject to project out regression matrix --
        #    (make errts like 3dDeconvolve, but more quickly)
        3dTproject -polort 0 -input pb04.$subj.r*.blursmall+tlrc.HEAD                  \
                   -mask mask_epi_anat_010.$subj+tlrc                              \
                   -censor censor_${subj}_combined_2.1D -cenmode ZERO             \
                   -dsort Local_FSWMe_rall+tlrc                                   \
                   -ort X_smoothed.nocensor.xmat.1D -prefix errts_smoothed.$subj.fanaticor

    # remove files to save disk_space
    rm *blursmall*

    end

end
