#!/bin/bash

FILE=$1

applywarp --in=$FILE.nii --out=${FILE}_fonovmni \
	  --ref=$HOME/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii \
	  --warp=$HOME/standard/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef --interp=spline
