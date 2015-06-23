#!/usr/local/bin/python

import csv
import corr

pop = open('SPECC_demographic.csv','r')
reader = csv.reader(pop)
pop_lines = [line for line in reader]
pop_header = pop_lines.pop(0)
#print pop_header
pop.close()

con = open('groupcov_ageexplore_n73.txt','r')
con_lines = [line.strip().split() for line in con.readlines()]
con_header = con_lines.pop(0)
#print con_header
con.close()

page = pop_header.index('AgeAtScan')
cage = con_header.index('age')
psex = pop_header.index('Sex')
csex = con_header.index('female')
plunaid = pop_header.index('LUNAID')
pmrid = pop_header.index('MRID')
clunaid = con_header.index('lunaid')

population = [{'id':i[pmrid],'age':float(i[page]),'sex':int(i[psex])-1} for i in pop_lines if i[plunaid] =="NA" and i[pmrid] not in corr.population_blacklist]
#print population

control_blacklist = [i[:5] for i in corr.control_blacklist]
#print control_blacklist
control = [{'id':i[clunaid],'age':float(i[cage]),'sex':int(i[csex])} for i in con_lines if i[clunaid] not in control_blacklist]
#print control


matches = []
for pt in population:
    #print pt
    potential_matches = [i for i in control if i['sex'] == pt['sex']]
    best = 0
    for i in range(1,len(potential_matches)):
        if abs(pt['age'] - potential_matches[i]['age']) < abs(pt['age'] - potential_matches[best]['age']):
            best = i
    matches.append(control.pop(control.index(potential_matches[best])))

print [i['id'] for i in matches]
