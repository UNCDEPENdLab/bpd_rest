#!/bin/bash

funcRuns=$(find $PWD -maxdepth 4 -iname brnswudktm_rest1_5.nii.gz)
set -ex
for f in $funcRuns; do
    cd $(dirname $f)

    #ROI_TempCorr.R -ts $f -rois /Volumes/Serena/SPECC/Neil/bpd_rest/neil/coordinate_generation/power264_mni2.3mm_appended.nii.gz -roi_reduce pca \
		   #-brainmask ~/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii \
		   #-censor $(dirname $f)/motion_info/censor_intersection.1D -out_file corr_rois_robust_new.txt -corr_method robust -fisherz

    ROI_TempCorr.R -ts $f -rois /Volumes/Serena/SPECC/Neil/bpd_rest/neil/coordinate_generation/power264_mni2.3mm_appended_culled.nii.gz -roi_reduce pca \
		   -brainmask ~/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii \
		   -censor $(dirname $f)/motion_info/censor_intersection.1D -out_file corr_rois_pearson_new_r_v2.txt -corr_method pearson

done
	 
