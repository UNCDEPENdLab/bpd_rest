#!/usr/local/bin/python
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import mni
import bct
from scipy import stats
import warnings


import corr

patient = "001RA_07DEC2013"
fileone = "corr_rois_pearson_new_r.txt"
#fileone = "corr_roimean_pearson.txt"

mat,mask = corr.load_patients([corr.population_folder+'/'+patient+'/'+fileone])
pat = mat[:,:,0]
#corr.pretty_print_2d(pat)
#print '=========='
hard_adj_mat = corr.map_adjacency_matrix(pat,corr.HARD,threshold=0.6)
hard_weighted_adj_mat = corr.map_adjacency_matrix(pat,corr.HARD_WEIGHTED,threshold=0.6)
hard_sparsed_adj_mat = corr.map_adjacency_matrix(pat,corr.HARD,percentile=0.9)
hard_weighted_sparsed_adj_mat = corr.map_adjacency_matrix(pat,corr.HARD_WEIGHTED,percentile=0.9)
soft_adj_mat = corr.map_adjacency_matrix(pat,corr.SOFT,beta=2)

adj_mat = [hard_adj_mat,hard_weighted_adj_mat,hard_sparsed_adj_mat,hard_weighted_sparsed_adj_mat,soft_adj_mat]
weighted = [0,1,0,1,1]
name = ["Hard threshold binary (equi-threshold)","Hard threshold weighted (equi-threshold)","Hard threshold binary (equi-sparse)","Hard threshold weighted (equi-sparse)","Soft / continuous weighted"]
#corr.pretty_print_2d(adj_mat[0])
#corr.draw_corr_matrix(adj_mat[0])

for i in range(0,len(adj_mat)):
    #print bct.bct.degrees_und(i)
    count = 0
    for j in adj_mat[i]:
        for k in j:
            if k != 0:
                count += 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Dangerous, but bct.py charpath is not well written and throws 0*inf warnings as written
        s = corr.network_measures(adj_mat[i],weighted=weighted[i])
    print name[i]
    print "# edges: %s"%(count/2)
    for k,v in s.iteritems():
        if type(v) ==np.ndarray:
            print k,"list of length %i" %len(v)
        else:
            print k,v
    #corr.draw_corr_matrix(adj_mat[i])
