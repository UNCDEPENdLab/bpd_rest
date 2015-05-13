import corr
import mni
from scipy import stats
import numpy as np
import bct
import warnings
import pylab

control_folder = corr.control_folder
control = corr.control[:10]
population_folder = corr.population_folder
population = corr.population[:10]
filename = 'corr_rois_pearson_new_r_v2.txt'

debug_timing = True
if debug_timing:
    import time
    currtime = time.clock()
names = [control,population]
everyone = [] # [0] is control supermatrix, [1] is population; format: [n][x][y][i]
everyone_mask = []
mat,mask = corr.load_patients([control_folder+'/'+i+'/'+filename for i in control])
everyone.append(mat)
everyone_mask.append(mask)
mat,mask = corr.load_patients([population_folder+'/'+i+'/'+filename for i in population])
everyone.append(mat)
everyone_mask.append(mask)
if debug_timing:
    newtime = time.clock()
    delta = newtime - currtime
    currtime = newtime
    print "load matrices",delta

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
if debug_timing:
    newtime = time.clock()
    delta = newtime - currtime
    currtime = newtime
    print "map matrices",delta

# generate stats
everyone_stats = [] # [n][dict]
for i in range(0,len(everyone_mapped)):
    """
    pop_stats = []
    for j in range(0,len(everyone_mapped[i])):
        pop_stats.append(corr.network_measures(everyone_mapped[i][j],weighted=weighted,limited=True))
    """
    func = corr.network_measures_helper_generator({'weighted':weighted,'limited':True})
    #pop_stats = map(func,everyone_mapped[i]) mapped version of above, serial
    parallel_func = corr.parallel_function(func)
    pop_stats = parallel_func(everyone_mapped[i],pool_size=2)

    everyone_stats.append(pop_stats)

if debug_timing:
    newtime = time.clock()
    delta = newtime - currtime
    currtime = newtime
    print "generate stats",delta
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

if debug_timing:
    newtime = time.clock()
    delta = newtime - currtime
    currtime = newtime
    print "calculate network measures",delta

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


if debug_timing:
    newtime = time.clock()
    delta = newtime - currtime
    currtime = newtime
    print "calculate measures",delta

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

if debug_timing:
    newtime = time.clock()
    delta = newtime - currtime
    currtime = newtime
    print "write matrices to file",delta
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
