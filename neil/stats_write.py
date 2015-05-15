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
filename = 'corr_rois_pearson_new_r_v2.txt'

output_folder = 'stats_output/'

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

# output format: one file per subject; each file is a 4D array [e][p][n][s] where e is edge-definition, p is a parameter [percentile, or beta], s is statistic, n is ROI

def write_subject(adj,techniques,filename):
    """
    Format in file: CSV with technique, parameter, statistic, roi
    input: 
        adjacency matrix (not yet mapped)
        techniques (array of tuples, including technique, percentile, and beta)
        filename to write to
    output: none (writes to file)
    """
    f = open(filename,'w')
    f.write('edge_def,parameter,roi,statistic,value\n')
    for t,p,b in techniques:
        mapped = corr.map_adjacency_matrix(adj,t,percentile=p,beta=b)
        weighted = False if t == corr.HARD else True
        stats = corr.network_measures(mapped,weighted=weighted,limited=True)
        local_measures = []
        for k,v in stats.iteritems():
            if type(v) == np.ndarray:
                local_measures.append(k)
        local_measures.sort()
        ET_map = {corr.HARD:"HARD",corr.HARD_WEIGHTED:'HARD_WEIGHTED',corr.SOFT:'SOFT'}
        parameter = p if p != -2 else b
        for i in range(0,len(stats[local_measures[0]])):
            for j in local_measures:
                f.write("{:s},{:f},{:d},{:s},{:f}\n".format(ET_map[t],parameter,i,j,stats[j][i]))
    f.close()

techniques = []
for i in [0.85,0.9,0.95,0.975,0.99]:
    techniques.append((corr.HARD,i,-2))
for i in [0.85,0.9,0.95,0.975,0.99]:
    techniques.append((corr.HARD_WEIGHTED,i,-2))
for i in [2,4,6,8,10,12,14,16]:
    techniques.append((corr.SOFT,-2,i))
for i in range(0,len(everyone)):
    for j in range(0,len(everyone[i][0,0,:])):
        write_subject(everyone[i][:,:,j],techniques,output_folder+names[i][j])
        if debug_timing:
            newtime = time.clock()
            delta = newtime - currtime
            currtime = newtime
            print "wrote file:",names[i][j],delta
