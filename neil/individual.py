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

#adj_mat = [hard_adj_mat,hard_weighted_adj_mat,hard_sparsed_adj_mat,hard_weighted_sparsed_adj_mat,soft_adj_mat]
adj_mat = [hard_adj_mat]
weighted = [0,1,0,1,1]
name = ["Hard threshold binary (equi-threshold)","Hard threshold weighted (equi-threshold)","Hard threshold binary (equi-sparse)","Hard threshold weighted (equi-sparse)","Soft / continuous weighted"]
short_name = ["hard_bin ET","hard_wei ET","hard_bin ES","hard_wei ES","power"]
#corr.pretty_print_2d(adj_mat[0])
#corr.draw_corr_matrix(adj_mat[0])
set_of_all_measures = set()
arr_measures = []
for i in range(0,len(adj_mat)):
    #print bct.bct.degrees_und(i)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Dangerous, but bct.py charpath is not well written and throws 0*inf warnings as written
        s = corr.network_measures(adj_mat[i],weighted=weighted[i])
    arr_measures.append(s)
    #for k,v in s.iteritems():
        #if k not in set_of_all_measures:
            #set_of_all_measures.add(k)
        #if type(v) ==np.ndarray:
            #print k,"list of length %i" %len(v)
        #else:
            #print k,v
    #corr.draw_corr_matrix(adj_mat[i])

corr.compare_networks(arr_measures,[short_name[0]])

#percentiles = [0.7,0.8,0.9,0.95,0.98]
#corr.compare_networks([corr.network_measures(corr.map_adjacency_matrix(pat,corr.HARD,percentile=i),0) for i in percentiles],[str(i) for i in percentiles])

patient = 0
mat = adj_mat[patient]
#mat,reorder_indices,cost = bct.bct.reorder_matrix(adj_mat[patient])
measures = arr_measures[patient]
print(measures['community_structure'])

graph,am = corr.create_graph(mat,weighted=False)
comm = graph.community_infomap().membership
print comm
poi_list = ["l_amygdala","r_amygdala","l_subgenual"]
read='ROI_nodes_new.node'
roi_list = corr.get_ROI_list(read)
poi_index = [mni.find_closest(mni.common_roi[i],roi_list)[0] for i in poi_list]
print poi_index
print [comm[i] for i in poi_index]

#colored_mat = [[comm[j] for j in range(0,len(i))] for i in mat]
#colored_mat = [[comm[j] if mat[i][j] != 0 else -2 for j in range(0,len(mat[i])) ] for i in range(0,len(mat))]
#corr.draw_corr_matrix(corr.sort2d(colored_mat,comm))

#write='my_ROI_new.node'
#corr.write_ROI_node_file(read,write,comm,poi_list)
