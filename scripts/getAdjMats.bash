#!/bin/bash
set -e

pipeline="mni_nosmooth_aroma"
fname="brnawuktm_rest1_5.nii.gz"

files=
for d in $( ls -d */ ); do
    if [ -f "${d}${pipeline}/rest1/${fname}" ]; then
	files="${files} $PWD/${d}${pipeline}/rest1/${fname}"
    fi
done

dirsrc=SPECC #or MMClock

#files=$( find $PWD -iname "brnswudktm_rest1_5.nii.gz" -ipath "*wavelet/rest1*" -type f | grep -v old )
#roifile=/gpfs/group/mnh5174/default/SPECC/Neil/bpd_rest/neil/coordinate_generation/power264_mni2.3mm_appended.nii.gz
roifile=/gpfs/group/mnh5174/default/lab_resources/parcellation/final_combined/7_Networks/Schaefer422_full_final.nii.gz
brainmask=/gpfs/group/mnh5174/default/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii
outdir=/gpfs/group/mnh5174/DEPENd_Box/Projects/RS_BPD_graph/adjmats_schaefer422_all
suffix=_schaefer422.txt.gz

set -x
for f in $files; do
    dir=$( dirname "$f" )
    id=$( echo $f | perl -pe "s:.*/${dirsrc}/MR_Proc/([^/]+)/mni.*:\1:" )
    while [ $( jobs -p | wc -l ) -ge 40 ]
    do
	sleep 20
    done
    
    if [ ! -f ${outdir}/${id}${suffix/_power269/_power269_pearson} ]; then

	ROI_TempCorr.R -ts "$f" -rois $roifile -roi_reduce pca -corr_method pearson,cor.shrink -pcorr_method pcor.shrink,ridge.net \
		       -brainmask $brainmask -njobs 1 -out_file ${outdir}/${id}${suffix} &

	#Since these are AROMA data, do not perform any scrubbing. Also, more recent ICA analyses do not suggest problems with steady state, so do not drop volumes
	#-censor $dir/motion_info/censor_intersection.1D -dropvols 5 \
	
	echo "" #space out log

    fi
done

wait
