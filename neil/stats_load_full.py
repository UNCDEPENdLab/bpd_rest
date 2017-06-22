#!/usr/local/bin/python
import math
import corr
import mni
import csv
import numpy as np
import pylab
from scipy import stats

folder = '/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_v2/'
controls = corr.control[:2]
population = corr.population[:2]
con = []
header_index = {}
header_type = {'edge_def':str,'parameter':float,'roi':int,'statistic':str,'value':float}
for i in controls:
    f = open(folder+i,'r')
    reader = csv.reader(f)
    header = reader.next()
    pat={}
    for j in range(0,len(header)):
        header_index[header[j]]=j
        pat[header[j]]=[]
    for line in reader:
        for j in range(0,len(line)):
            #row.append(header_type[header[j]](line[j]))
            pat[header[j]].append(header_type[header[j]](line[j]))
    f.close()
    con.append(pat)

print header_index
pat = con[0]
measures = []
edge_def = []
parameters = {}
#num_roi = np.max(pat['roi'])
roi = []
for i in range(0,len(pat['statistic'])):
    if pat['statistic'][i] not in measures:
        measures.append(pat['statistic'][i])
    if pat['edge_def'][i] not in edge_def:
        edge_def.append(pat['edge_def'][i])
        parameters[pat['edge_def'][i]]=[]
    if pat['roi'][i] not in roi:
        roi.append(pat['roi'][i])
    if pat['parameter'][i] not in parameters[pat['edge_def'][i]]:
        parameters[pat['edge_def'][i]].append(pat['parameter'][i])
    
print parameters

pylab.figure()
for k in range(0,len(measures)):
    for j in range(0,len(edge_def)):
        arr = [[] for i in range(0,len(roi))] # Rows: ROI, columns: statistic parameter
        for i in range(0,len(pat['statistic'])):
            if pat['statistic'][i] == measures[k] and pat['edge_def'][i]==edge_def[j]:
                arr[pat['roi'][i]].append(pat['value'][i])
        arr = np.array(arr)
        #pylab.subplot(len(measures),len(edge_def),k*len(edge_def)+j+1)
        ax = pylab.subplot(len(edge_def),len(measures),j*len(measures)+k+1)
        if j==0:
            pylab.title(measures[k])
        if k==0:
            pylab.ylabel(edge_def[j])
        pylab.pcolor(arr)
        pylab.colorbar()
        pylab.ylim([0,len(roi)])
        ax.set_xticks(np.arange(len(parameters[edge_def[j]]))+0.5, minor = False)
        ax.set_xticklabels(parameters[edge_def[j]],minor = False)
pylab.show()
# format of multidimensional array: [n][e][p][r][s]: n=patient,e=edge_def,p=parameter,r=roi,s=statistic
