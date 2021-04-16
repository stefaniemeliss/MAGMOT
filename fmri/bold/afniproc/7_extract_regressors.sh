#!/bin/bash

################################################################################
# non-linear registration to MNI template (non-linear 2009c)
################################################################################

# from https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/@SSwarper_sphx.html
# This script has dual purposes for processing a given subject's anatomical volume:
#    + to skull-strip the brain, and
#    + to calculate the warp to a reference template/standard space.
# Automatic snapshots of the registration are created, as well, to help the QC process.
# This program cordially ties in directly with afni_proc.py, so you can run it beforehand, 
# check the results, and then provide both the skull-stripped volume and the warps to the processing program.


# define path
DIR="/storage/shared/research/cinn/2018/MAGMOT"

# change directory to BIDS folder
BIDS_dir="$DIR"/MAGMOT_BIDS
# change directory to the raw NIFI files
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
subjects=($(ls -d sub-control*))
subjects=($(ls -d sub-experimental*))
#subjects=(sub-control001)

# define derivatives directory
deriv_dir="$DIR"/derivatives

# define which template to use and where to find them
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`

# detemine zeropad values
zp_60="sub-experimental010 sub-experimental014"
zp_70="sub-experimental004 sub-experimental005 sub-control007 sub-experimental008 sub-control011 sub-experimental016 sub-control021 sub-experimental024 sub-control027 sub-control039 sub-experimental040 sub-control043 sub-experimental048"
zp_80="sub-control001 sub-control002 sub-control003 sub-experimental006 sub-control013 sub-control015 sub-experimental018 sub-control019 sub-experimental022 sub-control025 sub-experimental028 sub-control029 sub-experimental030 sub-control031 sub-experimental032 sub-control033 sub-experimental034 sub-control035 sub-experimental036 sub-control037 sub-experimental044 sub-control045 sub-experimental046 sub-control049 sub-experimental050"
zp_90="sub-control009 sub-experimental012 sub-control017 sub-experimental020 sub-control023 sub-experimental026 sub-experimental038 sub-control041 sub-experimental042 sub-control047"

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

    echo "###################################################################################"	
    echo "@SSwarper for subject $subject"
    echo "###################################################################################"	

	# create output folder
	out_root=$deriv_dir/afniproc
    mkdir $out_root
    mkdir $out_root/$subject
	out_dir=$out_root/$subject/SSwarper
	mkdir $out_dir
    cd $out_dir

	# define BIDS anat folder
	anat_dir=$DIR/MAGMOT_BIDS/"${subject}"/anat
    anat_orig="$subject"_rec-NORM_T1w.nii.gz
    anat="$subject".nii.gz


    # copy anat
    3dcopy $anat_dir/$anat_orig ./$anat

    ########## REMOVE NECK AREA ##########

    # define prefix
    zp="$subject"_zp.nii.gz

    # define input for zero padding based on subject ID
	if echo $zp_60 | grep -w $subject > /dev/null; then 
       3dZeropad -I -60 -prefix $zp $anat
	elif echo $zp_70 | grep -w $subject > /dev/null; then 
       3dZeropad -I -70 -prefix $zp $anat
	elif echo $zp_80 | grep -w $subject > /dev/null; then 
       3dZeropad -I -80 -prefix $zp $anat
	elif echo $zp_90 | grep -w $subject > /dev/null; then 
       3dZeropad -I -90 -prefix $zp $anat
    fi

    ########## ALIGN CENTERS ##########

    # define prefix
    shft="$subject"_zp_shft.nii.gz

    # align centers
    @Align_Centers  \
        -cm \
        -prefix $shft     \
        -dset $zp \
        -base $template_path/$template

    ########## run @SSwarper ##########

	@SSwarper -input $shft  \
		-subid $subject     \
		-odir $out_dir      \
		-base $template_path/$template


    echo ""
    echo ""
    echo ""

    # as a last step, move anat_with_skull (anatUAC) anat anat_without_skull (anatSS) into derivatives folder

    # deine derivative dir
    anat_deriv=$DIR/derivatives/"${subject}"/anat
    #3dcopy anatUAC.$subject.nii $anat_deriv/"${subject}"_space-orig_res-anat_label-skulled_T1w.nii.gz
    3dcopy anatUAC.$subject.nii $anat_deriv/"${subject}"_label-withskull_T1w.nii.gz
    3dcopy anatSS.$subject.nii $anat_deriv/"${subject}"_label-skullstripped_T1w.nii.gz  



done

# this will go through the following steps (from @SSWarper -help)

  #1: Uniform-ize the input dataset's intensity via 3dUnifize.
  #     ==> anatU.sub007.nii
  # anatU.sub007.nii        = intensity uniform-ized original dataset (or, if '-unifize_off' used, a copy of orig dset);
  # anatUA.sub007.nii       = anisotropically smoothed version of the above (or, if '-aniso_off' used, a copy of anatU.*.nii)
  # anatUAC.sub007.nii      = ceiling-capped ver of the above (at 98%ile of non-zero values)
  #2: Strip the skull with 3dSkullStrip, with mildly agressive settings.
  #     ==> anatS.sub007.nii
  #3: Nonlinearly warp (3dQwarp) the result from #1 to the skull-on template, driving the warping to a medium level of refinement.
  #4: Use a slightly dilated brain mask from the template to crop off the non-brain tissue resulting from #3 (3dcalc).
  #5: Warp the output of #4 back to original anatomical space, along with the template brain mask, 
  #   and combine those with the output of #2 to get a better skull-stripped result in original space (3dNwarpApply and 3dcalc).
  #     ==> anatSS.sub007.nii
  #6  Restart the nonlinear warping, registering the output of #5 to the skull-off template brain volume (3dQwarp).
  #     ==> anatQQ.sub007.nii (et cetera)
  # anatQQ.sub007.nii       = skull-stripped dataset nonlinearly warped to the base template space;
  # anatQQ.sub007.aff12.1D  = affine matrix to transform original dataset to base template space;
  # anatQQ.sub007_WARP.nii  = incremental warp from affine transformation to nonlinearly aligned dataset;
  #7  Use @snapshot_volreg3 to make the pretty pictures.
  #     ==> AMsub007.jpg and MAsub007.jpg
  # AMsub007.jpg          = 3x3 snapshot image of the anatQQ.sub007.nii dataset with the edges from the base template overlaid -- to check the alignment;
  # MAsub007.jpg          = similar to the above, with the roles of the template and the anatomical datasets reversed.
  # QC_anatQQ.sub007.jp   = like AM*.jpg, but 3 rows of 8 slices
  # QC_anatSS.sub007.jpg  = check skullstripping in orig space: ulay is input dset, and olay is mask of skullstripped output (anatSS* dset)
  # init_qc_00_overlap_uinp_obase.jpg
  #  o [ulay] original source dset and [olay] original base dset
  #  o single image montage to check initial overlap of source and base, ignoring any obliquity that might be present (i.e., the way AFNIGUI does by default, and also how alignment starts)
  #  o if initial overlap is not strong, alignment can fail or produce weirdness
  #  o *if* either dset has obliquity, then an image of both after deobliquing with 3dWarp is created (*DEOB.jpg), and a text file about obliquity is also created (*DEOB.txt).








