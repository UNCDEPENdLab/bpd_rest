#!/bin/bash

#this volume contains mask numbers for each structure in the subcortical atlas
subcort=$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz

#lamygdala is 10, ramygdala is 20
3dcalc -a $subcort -expr 'amongst(a,10)' -prefix l_amygdala.nii.gz

#FSL MNI152 and Fonov 2009 MNI152 (official) are slightly different. Need to warp from FSL space (Harvard-Oxford) to
#Fonov 2009 space using warp coefficients computed previously.

# note: This will also transform the mask from 1mm voxels (given above maxprob map in 1mm vox) to
# 2.3mm voxels to match the resolution of the processed resting-state data.
applywarp --in=l_amygdala --out=l_amygdala_fonovmni \
	  --ref=$HOME/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii \
	  --warp=$HOME/standard/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef --interp=spline

#spline interpolation will generate fractional values for the mask.
#could use nearest neighbor interpolation (--interp=nn), but this tends to dilate the mask a bit more than is desirable
#FSL folks usually recommend thresholding the interpolated values at 0.5 to define the warped mask
fslmaths l_amygdala_fonovmni -thr 0.5 -bin l_amygdala_fonovmni


#same for r amygdala, and perhaps sgACC
3dclust -1Dformat -nosum -1dindex 0 -1tindex 0 -1thresh 0.5 -dxyz=1 1.01 10 l_amygdala_fonovmni.nii.gz > lamygdala_center.txt #CM LR CM PA CM IS

#use neurosynth for sgACC etc.

#location of FSL MNI template in case you want to overlay NeuroSynth ROIs onto structural
$HOME/standard/fsl_mni152/MNI152_T1_2mm.nii #MNI 152 T1 image in 2mm voxels


#append those coordinates to bb264 and number 265 etc.


#after all the appends, create the NIFTI mask volume
3dUndump -xyz -master /Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii -orient LPI -srad 5 \
	 -prefix power264_mni2.3mm_appended.nii.gz bb264coordinate_appended


