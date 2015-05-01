#!/usr/local/bin/python
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import copy
import pylab
import mni
import bct
from scipy import stats

base_folder = '/Volumes/Serena/Raj/Preprocess_Rest'

control_folder=base_folder
control=['10637_20140304', '10638_20140507', '10711_20140826', '10717_20140813', '10767_20140814', '10772_20140527', '10811_20140721', '10814_20140623', '10822_20140619', '10873_20140918', '10891_20140728', '10997_20140308', '11162_20140717', '11178_20140310', '11216_20141029', '11228_20140922', '11229_20140521', '11243_20140130', '11250_20140228', '11252_20140213', '11253_20140308', '11255_20140227', '11256_20140314', '11258_20140306', '11262_20140331', '11263_20140307', '11265_20141006', '11274_20140603', '11275_20140527', '11277_20140410', '11278_20140519', '11279_20140423', '11280_20140905', '11281_20140416', '11282_20141111', '11287_20140528', '11288_20140602', '11298_20140702', '11302_20140919', '11304_20140812', '11305_20140805', '11310_20140731', '11311_20140819', '11313_20140918', '11314_20140916', '11315_20140822', '11316_20140818', '11317_20140829', '11318_20140828', '11319_20140904', '11321_20140925', '11322_20140904', '11323_20141126', '11324_20141023', '11325_20141007', '11326_20140922', '11328_20141113', '11329_20141014', '11331_20141024', '11335_20141111', '11336_20141204', '11337_20141117', '11338_20141213', '11342_20150228', '11343_20141213', '11344_20141209', '11345_20141211', '11346_20150131', '11347_20141205', '11353_20150124'] # all controls
population_folder = base_folder+'/SPECC'
population= ['008JH_13JAN2014', '013jk_30Apr2014', '015cw_03May2014', '018LQ_26MAR2014', '019ec_04Aug2014', '020lr_03May2014', '023ds_07May2014', '025ay_10Jun2014', '027AD_18Sep2014', '031VN_09Sep2014', '037ll_25Aug2014', '038aa_03nov2014', '046ak_03Nov2014', '047ab_03nov2014', '048ah_18Dec2014', '050ai_06Nov2014', '054ls_12Jan2015', '058ab_15Jan2015', '059cr_08jan2015', '066dw_14Mar2015'] # removed first two (controls), last one is omitted as corr matrix had not yet completed
filename = 'corr_roimean_pearson.txt' # robust vs pearson

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
            else: #float
                print "{:5.3f}".format(arr[i][j]),
        print


"""
Loads each patient into the mat array (mat[x][y][i] where i is the patient, x and y represent correlation adjacency matrix for that patient (indexed by the Power et all 2011 264 node setup). Also creates a mask to indicate specific nodes where the correlation could not be calculated
Returns: matrix[][][],mask[][][]
"""
def load_patients(patient_list):
    mat = []
    mask = []
    for i in patient_list:
        pat = []
        patmask = []
        f = open(i,'r')
        for line in f:
            # TODO: Come up with suitable way to replace NA when the raw matrix is to be used in NBS testing
            l = [float(j.replace("NA","0")) for j in line.strip().split()] # 0 is a magic number signalling NA, needs to be dealt with when accessing the data later via the variable mask (set on next line); prev versions used arbitrary -42, but for NBS data processing 0 is a signal to exclude data from analysis so will be used
            ml = [True if j=="NA" else False for j in line.strip().split()]
            pat.append(l)
            patmask.append(ml)
        f.close()
        if len(mat) == 0:
            for x in range(0,len(pat)):
                mat.append([])
                mask.append([])
                for y in range(0,len(pat[x])):
                    mat[x].append([])
                    mask[x].append([])
        for x in range(0,len(pat)):
            for y in range(0,len(pat[x])):
                mat[x][y].append(pat[x][y])
                mask[x][y].append(patmask[x][y])
    return np.array(mat),np.array(mask)

"""
Creates an averaged 2D array over all patients, zeroes the diagonal. Requires
matrix[][][] as well as boolean mask[][][]. Returns averaged array[][].
Averages over all data points minus the masked out entries; will throw
exception if there is no valid element at a single node. Would caution using
with low number data sets as with even a few masked entries, the averaged data
point will be swung by as little as one remaining entry; may cause problems in
analysis later; todo: create a separate 2d array with total N for each data
point
"""
def average_patients(mat,mask):
    avg = []
    for x in range(0,len(mat)):
        avg.append([])
        for y in range(0,len(mat[x])):
            avg[x].append(0)
            count = 0
            for i in range(0,len(mat[x][y])):
                if mask[x][y][i]: # removes all items marked NA on original adjacency matrix
                    continue
                avg[x][y] += mat[x][y][i]
                count += 1
            avg[x][y]/=count # WILL throw exception if count = 0 due to all patients with masked data at single data point
            if x == y:
                avg[x][y] = 0 # Graph should not have loops
    return np.array(avg)



"""
Creates sparsified igraph Graph by putting all correlation data in a one dimensional
matrix and taking only highest percentile edges (tie density). Requires the
correlation matrix (recommended 0'd diagonal), percentile for cutoff (100 minus
tie-density), and takes optional argument to include weighted edges or binary.
Returns igraph Graph object and numpy array representing the sparsified matrix
"""
def create_graph(avg,percentile,weighted=True):
    G = ig.Graph()
    G.add_vertices(len(avg))
    sparsed = copy.deepcopy(avg)
    corr_list = []
    for i in range(0,len(avg)):
        for j in range(0,i):
            corr_list.append(avg[i][j])
    cl = np.array(corr_list)
    cutoff_percentile = 95 # percentile /  100 - tie density
    cutoff = np.percentile(cl,cutoff_percentile)
    #plt.hist(cl) # draws histogram of all correlations
    #plt.show()
    for i in range(0,len(avg)):
        for j in range(0,i):
            if avg[i][j] > cutoff:
                if weighted:
                    G.add_edge(i,j,weight=avg[i][j])
                else:
                    G.add_edge(i,j)
                    sparsed[i][j] = 1
                    sparsed[j][i] = 1
            else:
                sparsed[i][j] = 0
                sparsed[j][i] = 0
    return G,np.array(sparsed)

if __name__ == '__main__':
    #patients = population
    #folder = population_folder
    mat_control,mask_control = load_patients([control_folder+'/'+i+'/'+filename for i in control]) # will be mat[x][y][i] where i is trial, x and y represent the correlation adjacency matrix for that patient on the Power 264 node setup
    adj_control = average_patients(mat_control,mask_control)
    G_control,sparsed_adj_control = create_graph(adj_control,90,weighted=False)

    mat_population,mask_population = load_patients([population_folder+'/'+i+'/'+filename for i in population]) # will be mat[x][y][i] where i is trial, x and y represent the correlation adjacency matrix for that patient on the Power 264 node setup
    adj_population = average_patients(mat_population,mask_population) # adj_population[x][y] represents averaged adjacency matrix for all trials
    G_population,sparsed_adj_population = create_graph(adj_population,90,weighted=False)

    #ig.summary(G_control)
    #print G_control.degree(range(0,10))
    #print G_control.betweenness(range(0,10))
    #print G_control.edge_betweenness()[0:10]

    #draw_corr_matrix(adj_control)
    #draw_corr_matrix(adj_population)
    """
    Community detection algorithms
    can use optimal [too slow for 100+ node graphs, did not terminate after 45 min on 264 node], fastgreedy, infomap, and others
    """
    #mod = G_control.community_optimal_modularity() # TOO SLOW
    #membership = mod.membership
    #dendrogram = G_control.community_fastgreedy()
    #clusters = dendrogram.as_clustering()
    #membership = clusters.membership
    clusters_control = G_control.community_infomap()
    membership_control = clusters_control.membership
    #print clusters_control

    clusters_population = G_population.community_infomap()
    membership_population = clusters_population.membership
    #print clusters_population

    #print membership_control
    # ig.plot(cluster, vertex_label=range(0,len(adj_control)),vertex_label_size=8,bbox=[1000,1000]) # PLOT community clusters



    #Community comparators
    f1 = open('power_communities.txt','r')
    power = [int(f.strip()) for f in f1.read().split('\n') if len(f)>0]
    f1.close()
    print "NMI score (control to Power): ",ig.clustering.compare_communities(power,membership_control,method="nmi") # nmi, vi, etc
    print "NMI score (population to Power): ",ig.clustering.compare_communities(power,membership_population,method="nmi") # nmi, vi, etc
    print "NMI score (population to control): ",ig.clustering.compare_communities(membership_population,membership_control,method="nmi") # nmi, vi, etc


    f = open('ROI_nodes.node','r')
    roi = [line.strip().split('\t') for line in f]
    f.close()
    for i in range(0,len(membership_control)):
        roi[i][3] = membership_control[i]
    poi_index,dist = mni.find_closest(mni.common_roi['amygdala'],roi)

    
    print "PageRank (c): ",bct.bct.pagerank_centrality(adj_control,0.85)[poi_index]
    print "PageRank (p): ",bct.bct.pagerank_centrality(adj_population,0.85)[poi_index]
    #print "Weighted betweenness centrality (c): ",bct.bct.betweenness_wei(adj_control)[poi_index]
    #print "Weighted betweenness centrality (p): ",bct.bct.betweenness_wei(adj_population)[poi_index]
    #print "Binary betweenness centrality (c): ",bct.bct.betweenness_bin(sparsed_adj_control)[poi_index]
    #print "Binary betweenness centrality (p): ",bct.bct.betweenness_bin(sparsed_adj_population)[poi_index]
    print "Weighted clustering coefficient (c): ",bct.bct.clustering_coef_wu(adj_control)[poi_index]
    print "Weighted clustering coefficient (p): ",bct.bct.clustering_coef_wu(adj_population)[poi_index]
    #print "Binary clustering coefficient (c): ",bct.bct.clustering_coef_bu(sparsed_adj_control)[poi_index]
    #print "Binary clustering coefficient (p): ",bct.bct.clustering_coef_bu(sparsed_adj_population)[poi_index]

    #pretty_print_2d(mat_control[:,:,0])
    pagerank_control = [bct.bct.pagerank_centrality(mat_control[:,:,i],0.85)[poi_index] for i in range(0,len(mat_control[poi_index][poi_index]))]
    pagerank_population = [bct.bct.pagerank_centrality(mat_population[:,:,i],0.85)[poi_index] for i in range(0,len(mat_population[poi_index][poi_index]))]
    clustering_control = [bct.bct.clustering_coef_wu(mat_control[:,:,i])[poi_index] for i in range(0,len(mat_control[poi_index][poi_index]))]
    clustering_population = [bct.bct.clustering_coef_wu(mat_population[:,:,i])[poi_index] for i in range(0,len(mat_population[poi_index][poi_index]))]

    print stats.ttest_ind(pagerank_control,pagerank_population)
    print stats.ttest_ind(clustering_control,clustering_population)
    """
    #sorts the matrix by membership to more easily identify communities; in theory
    #the Power et al ROIs were selected and ordered such that large communities are
    #already sequential

    adj_control_sorted = sort2d(adj_control,membership_control)
    draw_corr_matrix(adj_control_sorted)
    """

    """
    membership2d = []
    for i in range(0,len(membership_control)):
        membership2d.append([])
        for j in range(0,len(membership_control)):
            if membership_control[i] == membership_control[j]:
                membership2d[i].append(membership_control[i])
            else:
                membership2d[i].append(-1)
    #draw_corr_matrix(membership2d)
    """

    """
    # write ROI node file with community membership data to be opened in BrainNet Viewer

    f = open('my_ROI.node','w')
    for i in range(0,len(roi)):
        for item in roi[i]:
            f.write(str(item)+'\t')
        f.write('\n')
    f.close()
    """
