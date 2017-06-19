#!/bin/bash

#this volume contains mask numbers for each structure in the subcortical atlas
subcort="$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz"

#fslview $subcort

#lamygdala is 10, ramygdala is 20 (roi_to_index.txt has my latest full list)
#3dcalc -a $subcort -expr 'amongst(a,10)' -prefix l_amygdala.nii.gz
#use subcortical_roi_to_nii.sh < roi_to_index.txt to generate the list automatically

#same for r amygdala, and perhaps sgACC
#3dclust -1Dformat -nosum -1dindex 0 -1tindex 0 -1thresh 0.5 -dxyz=1 1.01 10 l_amygdala.nii.gz > lamygdala_center.txt #CM LR CM PA CM IS

#use neurosynth for sgACC etc.

#location of FSL MNI template
fslview /Users/neil/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii /Volumes/Serena/Raj/Preprocess_Rest/power264_mni2.3mm.nii.gz & #MNI 152 T1 image in 2mm voxels
fslview $HOME/standard/fsl_mni152/MNI152_T1_2mm.nii #MNI 152 T1 image in 2mm voxels

#need to transform to official MNI 2009 here... will get code in a bit.

#append those coordinates to bb264 and number 265 etc.


#after all the appends, create the NIFTI mask volume
#3dUndump -xyz -master /Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii -orient LPI -srad 5 \
	 #-prefix power264_mni2.3mm_appended.nii.gz bb264coordinate_appended
