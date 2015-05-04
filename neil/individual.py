#!/usr/local/bin/python
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import mni
import bct
from scipy import stats

import corr

patient = "001RA_07DEC2013"
fileone = "corr_rois_pearson_new_r.txt"

mat,mask = corr.load_patients([corr.population_folder+'/'+patient+'/'+fileone])
pat = mat[:,:,0]
corr.pretty_print_2d(pat)

adj_mat = corr.map_adjacency_matrix(pat,corr.HARD_WEIGHTED,0.6)
corr.pretty_print_2d(adj_mat)

