#!/bin/bash
set -e

#[[ "$#" -lt 3 ]] && echo "Need 4 arguments: basedir, outdir, tempcorr_args, pipeline_suffix" && exit 1
basedir="/gpfs/group/mnh5174/default/MMClock/MR_Proc"
#basedir="/gpfs/group/mnh5174/default/SPECC/MR_Proc"
outdir="/gpfs/group/mnh5174/default/Michael/bpd_rest/roi_diagnostics_mnh"
tempcorr_args=""
pipeline_suffix="for_diagnostics"
njobs=20
dirsrc=$( basename $( dirname "$basedir" )) #should be MMClock or SPECC

cd "$basedir"

pipeline="mni_aroma_minimal_fsl"
fname="rnawuktm_rest1.nii.gz"

files=
for d in $( ls -d */ ); do
    if [ -f "${d}${pipeline}/rest1/${fname}" ]; then
	files="${files} $PWD/${d}${pipeline}/rest1/${fname}"
    fi
done

roifile=/gpfs/group/mnh5174/default/lab_resources/parcellation/final_combined_jul2018/Schaefer_422_final_jul2018.nii.gz #this is pre 95% thresholding

#brainmask=/gpfs/group/mnh5174/default/lab_resources/standard/fsl_mni152/MNI152_T1_2.3mm_brain_mask.nii
brainmask=/gpfs/group/mnh5174/default/Michael/bpd_rest/Schaefer422_95_groupmask.nii.gz

suffix=_schaefer422_fsl_95_groupmask.txt.gz
diag_suffix=_roidiagnostics_schaefer422_fsl_95_groupmask.csv

#set -x
for f in $files; do
    dir=$( dirname "$f" )
    id=$( echo $f | perl -pe "s:.*/${dirsrc}/MR_Proc/([^/]+)/mni.*:\1:" )
    while [ $( jobs -p | wc -l ) -ge ${njobs} ]
    do
	sleep 20
    done
    
    if [ ! -f ${outdir}/${id}${suffix/_schaefer422_fsl/_schaefer422_fsl_pearson} ]; then

        ROI_TempCorr.R -ts "$f" -rois $roifile -roi_reduce mean -corr_method pearson \
                       -brainmask $brainmask -njobs 1 -roi_diagnostics ${outdir}/${id}_${pipeline_suffix}${diag_suffix} \
                       -out_file ${outdir}/${id}_${pipeline_suffix}${suffix} ${tempcorr_args} > ${outdir}/${id}_tempcorr_${pipeline_suffix}.log 2>&1 &

    fi
done

wait
