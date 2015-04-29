#!/usr/local/bin/python
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import pylab

base_folder = '/Volumes/Serena/Raj/Preprocess_Rest'

control_folder=base_folder
control=['10637_20140304', '10638_20140507', '10711_20140826', '10717_20140813', '10767_20140814', '10772_20140527', '10811_20140721', '10814_20140623', '10822_20140619', '10873_20140918', '10891_20140728', '10997_20140308', '11162_20140717', '11178_20140310', '11216_20141029', '11228_20140922', '11229_20140521', '11243_20140130', '11250_20140228', '11252_20140213', '11253_20140308', '11255_20140227', '11256_20140314', '11258_20140306', '11262_20140331', '11263_20140307', '11265_20141006', '11274_20140603', '11275_20140527', '11277_20140410', '11278_20140519', '11279_20140423', '11280_20140905', '11281_20140416', '11282_20141111', '11287_20140528', '11288_20140602', '11298_20140702', '11302_20140919', '11304_20140812', '11305_20140805', '11310_20140731', '11311_20140819', '11313_20140918', '11314_20140916', '11315_20140822', '11316_20140818', '11317_20140829', '11318_20140828', '11319_20140904', '11321_20140925', '11322_20140904', '11323_20141126', '11324_20141023', '11325_20141007', '11326_20140922', '11328_20141113', '11329_20141014', '11331_20141024', '11335_20141111', '11336_20141204', '11337_20141117', '11338_20141213', '11342_20150228', '11343_20141213', '11344_20141209', '11345_20141211', '11346_20150131', '11347_20141205', '11353_20150124'] # all controls
population_folder = base_folder+'/SPECC'
population= ['001RA_07DEC2013', '005AI_06NOV2013', '008JH_13JAN2014', '013jk_30Apr2014', '015cw_03May2014', '018LQ_26MAR2014', '019ec_04Aug2014', '020lr_03May2014', '023ds_07May2014', '025ay_10Jun2014', '027AD_18Sep2014', '031VN_09Sep2014', '037ll_25Aug2014', '038aa_03nov2014', '046ak_03Nov2014', '047ab_03nov2014', '048ah_18Dec2014', '050ai_06Nov2014', '054ls_12Jan2015', '058ab_15Jan2015', '059cr_08jan2015', '066dw_14Mar2015']
filename = 'corr_roimean_pearson.txt' # robust vs pearson
patients = control
folder = control_folder

"""
Correlation matrix graphical generation
"""
def draw_corr_matrix(arr):
    arr_np = np.array(arr)
    pylab.pcolor(arr_np)
    pylab.colorbar()
    pylab.show()

"""
Returns the 2d array sorted by the index list
"""
def sort2d(arr,sorter_index):
    sorter = sorted(zip (sorter_index,range(0,len(sorter_index))))
    new_arr = []
    for i in range(0,len(arr)):
        new_arr.append([])
        for j in range(0,len(arr[i])):
            new_arr[i].append(arr[sorter[i][1]][sorter[j][1]])
    return new_arr

"""
Prints 2d array prettily, cutting off values as needed to fit in small space. Roughly assumes all rows have same number of columns
"""
def pretty_print_2d(arr,vsize = 10,hsize=10, integer = False):
    if len(arr) <= vsize:
        vertical_filler = False
    else:
        vertical_filler = True
    if len(arr[0]) <= hsize:
        horizontal_filler = False
    else:
        horizontal_filler = True
    for i in range(0,len(arr)):
        active_v_fill = False
        if vertical_filler and i >= vsize-2 and i<len(arr)-1:
            if i == vsize-2:
                active_v_fill = True
            else:
                continue
        for j  in range(0,len(arr[i])):
            if horizontal_filler and j >= hsize-2 and j < len(arr[i])-1:
                if j == hsize-2:
                    print " ....",
                continue
            if active_v_fill:
                print ".....",
                continue
            if integer:
                print "{:5d}".format(arr[i][j]),
            else:
                print "{:5.3f}".format(arr[i][j]),
        print

mat = [] # will be mat[i][x][y] where i is trial, x and y represent the correlation adjacency matrix for that patient on the Power 264 node setup
summed = [] # summed[x][y] represents averaged adjacency matrix for all trials

"""
Loads each patient into the mat array
"""
for i in patients:
    pat = []
    f = open(folder+'/'+i+'/'+filename,'r')
    for line in f:
        l = [float(j.replace("NA","0")) for j in line.strip().split()]
        pat.append(l)
    mat.append(pat)

#print len(mat), len(mat[0]), len(mat[0][0])

"""
Creates an averaged 2D array over all patients, zeroes the diagonal
"""
for x in range(0,len(mat[0])):
    summed.append([])
    for y in range(0,len(mat[0][x])):
        summed[x].append(0)
        for i in range(0,len(mat)):
            summed[x][y] += mat[i][x][y]
        summed[x][y]/=len(mat)
        if x == y:
            summed[x][y] = 0 # Graph should not have loops
        #print "{0:.2f} ".format(summed[x][y]),
    #print

draw_corr_matrix(summed)

"""
Creates igraph Graph by putting all correlation data in a one dimensional matrix and taking only highest percentile edges (tie density)
"""
G = ig.Graph()
G.add_vertices(len(summed))
corr_list = []
for i in range(0,len(summed)):
    for j in range(0,i):
        corr_list.append(summed[i][j])
cl = np.array(corr_list)
cutoff_percentile = 95 # percentile /  100 - tie density
cutoff = np.percentile(cl,cutoff_percentile)
#plt.hist(cl) # draws histogram of all correlations
#plt.show()
for i in range(0,len(summed)):
    for j in range(0,i):
        if summed[i][j] > cutoff:
            G.add_edge(i,j,weight=summed[i][j])
            #G.add_edge(i,j)

#ig.summary(G)
#print G.degree(range(0,10))
#print G.betweenness(range(0,10))
#print G.edge_betweenness()[0:10]

"""
Alternate adjacency matrix loading, loads whole graph (not pruned/thinned data)

G2 = ig.Graph.Read_Adjacency(folder+'/'+patients[0]+'/'+filename,attribute='weight',mode=ig.ADJ_UNDIRECTED) # full weighted graph
ig.summary(G2)
culled = G2.es.select(weight_gt = cutoff)["weight"]
"""

"""
Community detection algorithms
can use optimal [too slow for 100+ node graphs, did not terminate after 45 min on 264 node], fastgreedy, infomap, and others
"""
#mod = G.community_optimal_modularity() # TOO SLOW
#membership = mod.membership
#dendrogram = G.community_fastgreedy()
#clusters = dendrogram.as_clustering()
#membership = clusters.membership
clusters = G.community_infomap()
membership = clusters.membership
print clusters
#print membership
# ig.plot(cluster, vertex_label=range(0,len(summed)),vertex_label_size=8,bbox=[1000,1000]) # PLOT community clusters

membership2d = []
for i in range(0,len(membership)):
    membership2d.append([])
    for j in range(0,len(membership)):
        if membership[i] == membership[j]:
            membership2d[i].append(membership[i])
        else:
            membership2d[i].append(-1)
#draw_corr_matrix(membership2d)


#Community comparators
f1 = open('power_communities.txt','r')
power = [int(f.strip()) for f in f1.read().split('\n') if len(f)>0]
f1.close()
print ig.clustering.compare_communities(power,membership,method="nmi") # nmi, vi, etc


f = open('ROI_nodes.node','r')
roi = [line.strip().split('\t') for line in f]
f.close()

for i in range(0,len(membership)):
    roi[i][3] = membership[i]



# write ROI node file with community membership data to be opened in BrainNet Viewer
f = open('my_ROI.node','w')
for i in range(0,len(roi)):
    for item in roi[i]:
        f.write(str(item)+'\t')
    f.write('\n')
f.close()


#sorts the matrix by membership to more easily identify communities; in theory
#the Power et al ROIs were selected and ordered such that large communities are
#already sequential

#summed_sorted = sort2d(summed,membership)
#draw_corr_matrix(summed_sorted)
