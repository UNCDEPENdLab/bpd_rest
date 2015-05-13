import corr
import mni
from scipy import stats
import numpy as np
import bct
import warnings
import pylab

control_folder = corr.control_folder
control = corr.control
population_folder = corr.population_folder
population = corr.population
filename = 'corr_rois_pearson_new_r.txt'

names = [control,population]
everyone = [] # [0] is control supermatrix, [1] is population; format: [n][x][y][i]
everyone_mask = []
mat,mask = corr.load_patients([control_folder+'/'+i+'/'+filename for i in control])
everyone.append(mat)
everyone_mask.append(mask)
mat,mask = corr.load_patients([population_folder+'/'+i+'/'+filename for i in population])
everyone.append(mat)
everyone_mask.append(mask)

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
        pop.append(corr.map_adjacency_matrix(everyone[i][:,:,j],map_technique,percentile=percentile,beta=beta))
    everyone_mapped.append(pop)

# generate stats
everyone_stats = [] # [n][dict]
for i in range(0,len(everyone_mapped)):
    pop_stats = []
    for j in range(0,len(everyone_mapped[i])):
        pop_stats.append(corr.network_measures(everyone_mapped[i][j],weighted=weighted))
    everyone_stats.append(pop_stats)

global_measures = []
local_measures = []

# generate list of measures
for k,v in everyone_stats[0][0].iteritems():
    if type(v) == np.ndarray:
        local_measures.append(k)
    else:
        global_measures.append(k)

global_measures.remove("num_edges")
global_measures.remove("num_vertices")
global_measures.remove("mean_degree")
global_measures.remove("assortativity_binary")
global_measures.remove("mean_assortativity_binary")
local_measures.remove("eccentricity")
local_measures.remove("community_structure")
if map_technique==corr.SOFT:
    global_measures.remove("mean_assortativity_weighted")
    global_measures.remove("assortativity_weighted")
    global_measures.remove("global_efficiency")
    global_measures.remove("giant_component")
    global_measures.remove("diameter")
    global_measures.remove("radius")
    global_measures.remove("charpath")
    local_measures.remove("degree")
    local_measures.remove("betweenness_binary")

# remainder of code assumes that everyone.* arrays have only 2 entries, control = 0, pop = 1

# compare global measures
print "{:30s} {:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}".format("key","mean_c","mean_p","std_c","std_p","t-score","p-value")
for key in global_measures:
    control = [i[key] for i in everyone_stats[0]]
    population = [i[key] for i in everyone_stats[1]]
    s = stats.ttest_ind(control,population)
    print "{:30s} {:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}".format(key[:30], np.mean(control),np.mean(population),np.std(control),np.std(population),s[0],s[1])

# compare local measures at specified ROIs
rois = ["l_amygdala","r_amygdala","l_subgenual"]
read='ROI_nodes_new_v2.node'
roi_list = corr.get_ROI_list(read)
roi_index = [mni.find_closest(mni.common_roi[i],roi_list)[0] for i in rois]

for i in range(0,len(roi_index)):
    roi = rois[i]
    ind = roi_index[i]
    print ""
    print "ROI: ",roi,", node",ind
    print "{:30s} {:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}".format("key","mean_c","mean_p","std_c","std_p","t-score","p-value")
    for key in local_measures:
        control = [j[key][ind] for j in everyone_stats[0]]
        population = [j[key][ind] for j in everyone_stats[1]]
        s = stats.ttest_ind(control,population)
        print "{:30s} {:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}".format(key[:30], np.mean(control),np.mean(population),np.std(control),np.std(population),s[0],s[1])



"""
# write out all local measures in 1 file per subject, 1 line per ROI
all_measures = local_measures
out_folder = '/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_bin/'
# create file for each subject, rows are ROIs, columns are statistical measures
for i in range(0,len(everyone_stats)):
    for j in range(0,len(everyone_stats[i])):
        out_filename = names[i][j]
        f = open(out_folder+out_filename,'w')
        f.write("ROI,"+",".join(all_measures)+'\n')
        for k in range(0,len(everyone_stats[i][j][all_measures[0]])):
            f.write(str(k)+','+','.join([str(everyone_stats[i][j][measure][k]) for measure in all_measures])+'\n')
        print out_filename
        f.close()
"""

"""
roi = 0
measures = ["strength","local_efficiency","clustering_weighted","pagerank","betweenness_binary","eigenvector_centrality_und"]
for i in range(0,len(measures)):
    pylab.subplot(2,3,i+1)
    pylab.title(measures[i])
    pylab.hist([j[measures[i]][roi_index[roi]] for j in everyone_stats[0]],20)
pylab.show()
"""
