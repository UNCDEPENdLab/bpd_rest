#!/usr/local/bin/python
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import mni
import bct
from scipy import stats

import corr

patient = "001RA_07DEC2013"
fileone = "corr_rois_pearson_new.txt"

mat,mask = corr.load_patients([corr.population_folder+'/'+patient+'/'+fileone])
corr.pretty_print_2d(mat[:,:,0])

pat = mat[:,:,0]
