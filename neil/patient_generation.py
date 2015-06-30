#!/usr/local/bin/python
"""
Finds patients who show movement above threshold for a sufficient amount of time of the run
"""
import glob

base_folder = '/Volumes/Serena/Raj/Preprocess_Rest'
control_folder = ''
patient_folder = 'SPECC'
fd_file = 'motion_info/fd.txt'
threshold = 0.5
cutoff_ratio = 0.15

controls = glob.glob(base_folder+'/'+control_folder+'/*/'+fd_file)
patients = glob.glob(base_folder+'/'+patient_folder+'/*/'+fd_file)

all_pat = [controls,patients]

results = {}
for pop in all_pat:
    for pat in pop:
        f = open(pat,'r')
        lines = [float(line.strip()) for line in f.readlines()]
        f.close()
        s = sum([True if i>threshold else False for i in lines])
        ratio = float(s)/len(lines)
        results[pat] = ratio
        print "{:<80}\t{}\t{:2.3f} {}".format(pat,s,ratio,True if ratio > cutoff_ratio else "")

