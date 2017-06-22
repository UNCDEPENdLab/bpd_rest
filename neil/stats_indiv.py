#!/usr/local/bin/python
import math
import corr
import csv
import numpy as np
import pylab

folder = '/Volumes/Serena/Raj/Preprocess_Rest/SPECC/'
patient = '008JH_13JAN2014'
files = '/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_bin/'+patient
filename = 'corr_rois_pearson_new_r_v2.txt'
datafile = folder+patient+'/'+filename

f = open(files,'rt')
reader = csv.reader(f)
data = [r for r in reader]
measures = data.pop(0)
size = len(measures)
f.close()

#pat,mask = corr.load_patients([datafile])
#pat = pat[:,:,0]

data = np.array([[float(j) for j in i] for i in data])

for i in range(1,size): # skip ROI column
    for j in range(1,size):
        pylab.subplot(size-1,size-1,(i-1)*(size-1)+j)
        if i==1:
            pylab.title(measures[j])
        if j==1:
            pylab.ylabel(measures[i])
        pylab.scatter(data[:,j],data[:,i])
pylab.show()

pylab.plot()
for i in range(1,size):
    pylab.subplot(int(math.ceil((size-1)/math.ceil(math.sqrt((size-1))))),int(math.ceil(math.sqrt((size-1)))),i)
    pylab.hist(data[:,i],50)
    pylab.title(measures[i])
pylab.show()
