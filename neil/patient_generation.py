#!/usr/local/bin/python
"""
Finds patients who show movement above threshold for a sufficient amount of time of the run
"""
import glob
import os

fd_file = 'motion_info/fd.txt'
allsubj = open('../subject_list.txt','r')
subjs = [os.path.split(i.strip())[0]+'/'+fd_file for i in allsubj.readlines()]
allsubj.close()
threshold = 0.5
cutoff_ratio = 0.20

results = {}
for pat in subjs:
    f = open(pat,'r')
    lines = [float(line.strip()) for line in f.readlines()]
    f.close()
    s = sum([True if i>threshold else False for i in lines])
    ratio = float(s)/len(lines)
    results[pat] = ratio
    print "{:<80}\t{}\t{:2.3f} {}".format(pat,s,ratio,True if ratio > cutoff_ratio else "")

