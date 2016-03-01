#!/usr/local/bin/python
import corr
import bct
import os
import numpy as np
from scipy import stats
import sys

#control_folder = corr.control_folder
#control_folder_post = ''
#population_folder = corr.population_folder
#population_folder = '/Volumes/Serena/SPECC/MR_Proc_Rest/'
#population_folder_post = 'mni_5mm_wavelet/rest1/'
allsubj = open('../subject_list.txt','r')
subjs = [i.strip() for i in allsubj.readlines()]
allsubj.close()
filename = 'corr_rois_pearson_new_r_v2.txt'
#blacklist = ["001RA_07DEC2013", "005AI_06NOV2013", "023ds_07May2014", "050ai_06Nov2014", "0531lw_16Dec2014"] # replaced by blacklist.txt file

matched_controls_file = open('demographics/limited_controls.txt','r')
matched_controls = [i.strip() for i in matched_controls_file.readlines()]
matched_controls_file.close()

blacklist_file = open('../blacklist.txt','r')
blacklist = [i.strip().split()[0] for i in blacklist_file.readlines()]
blacklist_file.close()

roi_file = open('coordinate_generation/bb264coordinate_appended_culled_labels.txt','r')
rois = [r.strip('"').split('\t') for r in roi_file.readlines()][1:]

all_controls = [i for i in subjs if i.find("SPECC")==-1]
all_population = [i for i in subjs if i.find("SPECC")!=-1]

controls = [i for i in all_controls if i not in blacklist]
#control = matched_controls
population = [i for i in all_population if i not in blacklist]

everyone = [] # [0] is control supermatrix, [1] is population; format: [n][x][y][i]
everyone_mask = []
mat,mask = corr.load_patients(controls)
everyone.append(mat)
everyone_mask.append(mask)
mat,mask = corr.load_patients(population)
everyone.append(mat)
everyone_mask.append(mask)

print 'Significant edges:'
for i in range(0,len(everyone[0])):
    for j in range(0,i):
        con = [np.arctanh(k) for k in everyone[0][i][j]]
        pop = [np.arctanh(k) for k in everyone[1][i][j]]
        s = stats.ttest_ind(con,pop)
        if abs(s[0]) > 4.3:
            print i,j,s[0],s[1],rois[i][1],rois[j][1]


# mapping the adjacency matrices
everyone_mapped = [] # WARNING: changes format from everyone; [n][i][x][y] where n is population as above, i is pt, x,y are adjacency matrix
percentile,beta = -2,-2
# comment out the undesired methods
map_technique,percentile = corr.HARD,0.9
#map_technique,percentile = corr.HARD_WEIGHTED,0.9
#map_technique,beta = corr.SOFT,12
weighted = False if map_technique == corr.HARD else True
for i in range(0,len(everyone)): # pick group (control vs pop)
    pop = []
    for j in range(0,len(everyone[i][0,0,:])): # pick patient
        pop.append(corr.map_adjacency_matrix(everyone[i][:,:,j],map_technique,percentile=percentile,beta=beta,ignore_negative_weights=True))
    everyone_mapped.append(pop)

# generate stats
everyone_stats = [] # [n][dict]
for i in range(0,len(everyone_mapped)):
    func = corr.network_measures_helper_generator({'weighted':weighted,'limited':True}) # limited limits data collection to fewer measures, particularly filtering ones out that take a long time and don't seem to have good function for this application
    pop_stats = map(func,everyone_mapped[i]) # mapped version of prior for loop, serial, works well
    #parallel_func = corr.parallel_function(func) # DOES NOT WORK CURRENTLY, likely due to memory inefficiency in the centrality measure calculations in BCT
    #pop_stats = parallel_func(everyone_mapped[i],pool_size=2)

    everyone_stats.append(pop_stats)

global_measures = []
local_measures = []

# generate list of measures
for k,v in everyone_stats[0][0].iteritems():
    if type(v) == np.ndarray:
        local_measures.append(k)
    else:
        global_measures.append(k)

if "num_edges" in global_measures: global_measures.remove("num_edges")
if "num_vertices" in global_measures: global_measures.remove("num_vertices")
if "mean_degree" in global_measures: global_measures.remove("mean_degree")
if "assortativity_binary" in global_measures: global_measures.remove("assortativity_binary")
if "mean_assortativity_binary" in global_measures: global_measures.remove("mean_assortativity_binary")
if "eccentricity" in local_measures: local_measures.remove("eccentricity")
if "community_structure" in local_measures: local_measures.remove("community_structure")
if map_technique==corr.SOFT:
    if "mean_assortativity_weighted" in     global_measures:     global_measures.remove("mean_assortativity_weighted")
    if "assortativity_weighted" in global_measures:     global_measures.remove("assortativity_weighted")
    if "global_efficiency" in global_measures:     global_measures.remove("global_efficiency")
    if "giant_component" in global_measures:     global_measures.remove("giant_component")
    if "diameter" in global_measures:     global_measures.remove("diameter")
    if "radius" in global_measures:     global_measures.remove("radius")
    if "charpath" in global_measures:     global_measures.remove("charpath")
    if "degree" in local_measures:     local_measures.remove("degree")
    if "betweenness_binary" in local_measures:     local_measures.remove("betweenness_binary")
local_measures.sort()
global_measures.sort()


roi_index = range(0,len(everyone[0]))
#roi_index = [30, 190,212,228,230,234]

print ""
print "Global measures:"
for key in global_measures:
    con = [j[key] for j in everyone_stats[0]]
    pop = [j[key] for j in everyone_stats[1]]
    s = stats.ttest_ind(con,pop)
    print "{:30s} {:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}".format(key[:30],np.mean(con),np.mean(pop),np.std(con),np.std(pop),s[0],s[1])

for ind in roi_index:
    printstr = ""
    do_i_print = False
    printstr+= "node "+str(ind)+" "+rois[ind][1]+'\n'
    printstr+="{:30s} {:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}\n".format("key","mean_c","mean_p","std_c","std_p","t-score","p-value")
    for key in local_measures:
        con = [j[key][ind] for j in everyone_stats[0]]
        pop = [j[key][ind] for j in everyone_stats[1]]
        s = stats.ttest_ind(con,pop)
        if s[1] < 0.002:
            do_i_print = True
            printstr+="{:30s} {:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}\n".format(key[:30], np.mean(con),np.mean(pop),np.std(con),np.std(pop),s[0],s[1])
    if do_i_print:
        print ""
        print printstr
