import corr
import mni
from scipy import stats
import numpy as np
import bct
import warnings

control_folder = corr.control_folder
control = corr.control
population_folder = corr.population_folder
population = corr.population
filename = 'corr_rois_pearson_new_r.txt'

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
everyone_stats = []
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
local_measures.remove("eccentricity")

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
read='ROI_nodes_new.node'
roi_list = corr.get_ROI_list(read)
roi_index = [mni.find_closest(mni.common_roi[i],roi_list)[0] for i in rois]

for i in range(0,len(rois)):
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
