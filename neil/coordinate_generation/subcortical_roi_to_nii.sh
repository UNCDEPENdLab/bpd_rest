#!/bin/bash

#this volume contains mask numbers for each structure in the subcortical atlas
subcort="$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz"
output=subcortical_MNI.txt #only gets addended, not overwritten

while read -r line
do
	name=$(echo $line | awk '{print $1;}')
	index=$(echo $line | awk '{print $2;}')
	3dcalc -a $subcort -expr "amongst(a,$index)" -prefix $name.nii.gz
	printf "%s " "$name" >> $output
	3dclust -1Dformat -nosum -1dindex 0 -1tindex 0 -1thresh 0.5 -dxyz=1 1.01 10 $name.nii.gz | tail -n1 | awk '{print $2,$3,$4;}'  >> $output
	rm $name.nii.gz
done
#3dcalc -a $subcort -expr 'amongst(a,10)' -prefix l_amygdala.nii.gz
