#!/usr/local/bin/python
import glob

base_folder = '/Volumes/Serena/Raj/Preprocess_Rest'
control_folder = ''
patient_folder = 'SPECC'
fd_file = 'motion_info/fd.txt'

controls = glob.glob(base_folder+'/'+control_folder+'/*/'+fd_file)
patients = glob.glob(base_folder+'/'+patient_folder+'/*/'+fd_file)

all_pat = [controls,patients]

results = {}
for pop in all_pat:
    for pat in pop:
        f = open(pat,'r')
        lines = [float(line.strip()) for line in f.readlines()]
        f.close()
        threshold = 0.5
        s = sum([True if i>threshold else False for i in lines])
        ratio = float(s)/len(lines)
        results[pat] = ratio
        print "{:<80}\t{}\t{:2.3f} {}".format(pat,s,ratio,True if ratio > 0.15 else "")

