#!/usr/local/bin/python
import math
import corr
import mni
import csv
import numpy as np
import pylab
from scipy import stats
from IPython import embed

def permutation_test(pooled,sizeZ,sizeY,one_sided = False):
    np.random.shuffle(pooled)
    starZ = pooled[:sizeZ]
    starY = pooled[-sizeY:]
    return starZ.mean() - starY.mean() if one_sided else abs(starZ.mean() - starY.mean())

"""
Runs Monte Carlo simulation of the permutation test
Input:
    Control np.array
    Population np.array
    num_samples (default 10000)
    one_sided (default False)
Output:
    percentile of values outside of the given difference of means (if two-sided, will return values w/ absolute value greater than true diff
"""
def run_permutation_monte_carlo(control,population,num_samples=10000,one_sided=False,hist=False):
    pooled = np.hstack([control,population])
    delta = control.mean() - population.mean()
    estimates = np.array(map(lambda x: permutation_test(pooled,control.size,population.size,one_sided),range(num_samples)))
    diff_count = len(np.where(estimates <=delta)[0]) 
    hat_asl_perm = 1.0 - (float(diff_count)/float(num_samples))
    if hist:
        h,b,patch = pylab.hist(estimates,50,normed=True)
        patch[np.where(b > delta)[0][0]].set_color('r')
    return hat_asl_perm



folder = '/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_bin/'
controls = corr.control
population = corr.population
con = []
pop = []
measures = []
for i in controls:
    f = open(folder+i,'r')
    reader = csv.reader(f)
    data = [r for r in reader]
    measures = data.pop(0)
    data = np.array([[float(j) for j in i] for i in data])
    f.close()
    con.append(data)
for i in population:
    f = open(folder+i,'r')
    reader = csv.reader(f)
    data = [r for r in reader]
    measures = data.pop(0)
    data = np.array([[float(j) for j in i] for i in data])
    f.close()
    pop.append(data)
con = np.array(con)
pop = np.array(pop)
rois = ["l_amygdala","r_amygdala","l_subgenual"]
read='ROI_nodes_new_v2.node'
roi_list = corr.get_ROI_list(read)
roi_index = [mni.find_closest(mni.common_roi[i],roi_list)[0] for i in rois]

"""for i in range(0,len(roi_index)):
    roi = rois[i]
    ind = roi_index[i]
    print ""
    print "ROI: ",roi,", node",ind
    print "{:30s} {:>10s}{:>10s}{:>10s}{:>10s}{:>10s}{:>10s}".format("key","mean_c","mean_p","std_c","std_p","t-score","p-value")
    for x in range(1,len(measures)):
        key = measures[x]
        c = con[:,ind,x]
        p = pop[:,ind,x]
        s = stats.ttest_ind(c,p)
        print "{:30s} {:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}".format(key[:30], np.mean(c),np.mean(p),np.std(c),np.std(p),s[0],s[1])
"""

pylab.figure()
for i in range(0,len(roi_index)):
    roi = rois[i]
    ind = roi_index[i]
    print roi
    for j in range(1,len(measures)):
        key = measures[j]
        c = con[:,ind,j]
        p = pop[:,ind,j]
        pylab.subplot(len(roi_index),len(measures)-1,i*(len(measures)-1)+j)
        pylab.title(key)
        if j == 1:
            pylab.ylabel(roi)
        monte = run_permutation_monte_carlo(c,p,one_sided=True,hist=True)
        print "{:30s} {:10.4f}".format(key[:30],monte)

pylab.show()
embed()
