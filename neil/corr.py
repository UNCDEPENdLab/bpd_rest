#!/usr/local/bin/python
# corr.py: Process and analyze correlation/adjacency matrices generated from resting state fMRI
# @author: Neil Munjal
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import copy
import pylab
import mni
import bct
from scipy import stats
from functools import partial

base_folder = '/Volumes/Serena/Raj/Preprocess_Rest'

control_folder=base_folder
control=['10637_20140304', '10638_20140507', '10711_20140826', '10717_20140813', '10767_20140814', '10772_20140527', '10811_20140721', '10814_20140623', '10822_20140619', '10873_20140918', '10891_20140728', '10997_20140308', '11162_20140717', '11178_20140310', '11216_20141029', '11228_20140922', '11229_20140521', '11243_20140130', '11250_20140228', '11252_20140213', '11253_20140308', '11255_20140227', '11256_20140314', '11258_20140306', '11262_20140331', '11263_20140307', '11265_20141006', '11274_20140603', '11275_20140527', '11278_20140519', '11279_20140423', '11280_20140905', '11281_20140416', '11282_20141111', '11287_20140528', '11288_20140602', '11298_20140702', '11302_20140919', '11304_20140812', '11305_20140805', '11310_20140731', '11311_20140819', '11313_20140918', '11314_20140916', '11315_20140822', '11316_20140818', '11317_20140829', '11318_20140828', '11319_20140904', '11321_20140925', '11322_20140904', '11323_20141126', '11324_20141023', '11325_20141007', '11326_20140922', '11328_20141113', '11329_20141014', '11331_20141024', '11335_20141111', '11336_20141204', '11337_20141117', '11338_20141213', '11342_20150228', '11343_20141213', '11344_20141209', '11345_20141211', '11346_20150131', '11347_20141205', '11353_20150124'] # all controls except '11277_20140410' because it was missing at last generation
control_blacklist = ['10637_20140304','10711_20140826','11305_20140805'] # Too much motion artifact as noted by patient_generation.py

population_folder = base_folder+'/SPECC'
population= ['008JH_13JAN2014', '013jk_30Apr2014', '015cw_03May2014', '018LQ_26MAR2014', '019ec_04Aug2014', '020lr_03May2014', '023ds_07May2014', '025ay_10Jun2014', '027AD_18Sep2014', '031VN_09Sep2014', '037ll_25Aug2014', '038aa_03nov2014', '046ak_03Nov2014', '047ab_03nov2014', '048ah_18Dec2014', '049tm_17Apr2015', '050ai_06Nov2014', '0531lw_16Dec2014', '054ls_12Jan2015', '057as-09Dec2014', '058ab_15Jan2015', '059cr_08jan2015', '066dw_14Mar2015', '067sm_23Apr2015', '071eh_09Apr2015'] # removed first two (controls); TODO: dynamically generate this (and control) list and add the first two items in folder to blacklist
population_blacklist = ['023ds_07May2014','050ai_06Nov2014','0531lw_16Dec2014','028th_07Jul2014'] # First three: too much motion artifact as noted by patient_generation.py; 4th: rest data doesn't exist
#filename = 'corr_roimean_pearson.txt' # robust vs pearson
filename = 'corr_rois_pearson_new_r_v2.txt' # robust vs pearson; r = using r values, instead of z/t-scores as was previously done; v2 due to minor adjustments made to preprocessing script

def filter_list(orig,blacklist):
    """
    Returns the list with the blacklisted items deleted; ignores items not found in orig
    """
    return [i for i in orig if i not in blacklist]

def draw_corr_matrix(arr,show=True):
    """
    Correlation matrix graphical generation
    """
    arr_np = np.array(arr)
    pylab.pcolor(arr_np)
    pylab.colorbar()
    if show:
        pylab.show()

def sort2d(arr,sorter_index):
    """
    Returns the 2d array sorted by the index list; primary use is to sort communities together
    """
    sorter = sorted(zip (sorter_index,range(0,len(sorter_index))))
    new_arr = []
    for i in range(0,len(arr)):
        new_arr.append([])
        for j in range(0,len(arr[i])):
            new_arr[i].append(arr[sorter[i][1]][sorter[j][1]])
    return new_arr

def pretty_print_2d(arr,vsize = 10,hsize=10, integer = False):
    """
    Prints 2d array prettily, cutting off values as needed to fit in small space. Roughly assumes all rows have same number of columns
    """
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
        if not active_v_fill:
            print "{:4d}: ".format(i),
        else:
            print "....: ",
        for j  in range(0,len(arr[i])):
            if horizontal_filler and j >= hsize-2 and j < len(arr[i])-1:
                if j == hsize-2:
                    print " .....",
                continue
            if active_v_fill:
                print "......",
                continue
            if integer:
                if arr[i][j] == 0:
                    print "".ljust(6),
                else:
                    print "{:6d}".format(arr[i][j]),
            else: #float
                if arr[i][j] == 0:
                    print "{:6d}".format(0),
                else:
                    print "{: 5.3f}".format(arr[i][j]),
        print


def load_patients(patient_list):
    """
    Loads each patient into the mat array (mat[x][y][i] where i is the patient, x and y represent correlation adjacency matrix for that patient (indexed by the Power et all 2011 264 node setup). Also creates a mask to indicate specific nodes where the correlation could not be calculated
    Returns: matrix[][][],mask[][][]
    """
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

def average_patients(mat,mask):
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



def create_graph(avg,percentile=0,weighted=True):
    """
    Creates sparsified igraph Graph by putting all correlation data in a one dimensional
    matrix and taking only highest percentile edges (tie density). Requires the
    correlation matrix (recommended 0'd diagonal), percentile for cutoff (100 minus
    tie-density), and takes optional argument to include weighted edges or binary.
    Returns igraph Graph object and numpy array representing the sparsified matrix

    Note: can be used to generate a graph with all edges intact if percentile is given as 0
    """
    G = ig.Graph()
    G.add_vertices(len(avg))
    sparsed = copy.deepcopy(avg)
    corr_list = []
    for i in range(0,len(avg)):
        for j in range(0,i):
            corr_list.append(avg[i][j])
    cl = np.array(corr_list)
    cutoff = np.percentile(cl,percentile)
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

HARD = 0
HARD_WEIGHTED = 1
SOFT = 3
def map_adjacency_matrix(mat, map_type, threshold = -2, beta = -2, percentile = -2):
    """
    Creates a mapped adjacency matrix.

    Input: 
        mat: 2D n x n adjacency matrix, preferably n range [-1,1] (only important for
            continuous power-law distribution
        map_type: HARD (default)| HARD_WEIGHTED | SOFT
            HARD applies a hard threshold (above which an unweighted edge is placed)
            HARD_WEIGHTED applies threshold (above which a weighted edge has weight = r)
            SOFT applies continuous power function ((r+1)/2)^beta to map r from [-1,1] to [0,1]
                beta of  1 degenerates to simple linear weighting (but mapping -1,1 to 0,1)
        threshold (required if using HARD or HARD_WEIGHTED)
        beta (required if using SOFT)
        percentile: if defined, will supercede threshold; range [0-1]. Will return
            a hard/hard_weighted matrix with approximately ((1-percentile)*100)% of
            original edges; 1-percentile is tie-density. Can be used to generate a
            "equi-sparse" networks across patients, as opposed to "equi-threshold"

    Output:
        adj: 2D n x n adjacency matrix; if HARD, it is binary, else it is float with range [0,1] (only guranteed if input data is [-1,1]
    """
    global HARD, HARD_WEIGHTED, SOFT
    if ( (map_type == HARD or map_type == HARD_WEIGHTED) and threshold == -2 and percentile == -2):
        raise TypeError("Hard thresholded graphs require either a threshold or a percentile")
    if (map_type == SOFT and beta < 1):
        raise TypeError("Soft thresholded graphs require a beta from [1,:]")
    adj = copy.deepcopy(mat)

    if (map_type == HARD or map_type == HARD_WEIGHTED):
        cutoff = threshold
        if percentile != -2:
            corr_list = []
            for i in range(0,len(mat)):
                for j in range(0,i):
                    corr_list.append(mat[i][j])
            cl = np.array(corr_list)
            cutoff_percentile = percentile * 100
            cutoff = np.percentile(cl,cutoff_percentile)
        #print cutoff
        for i in range(0,len(mat)):
            for j in range(0,len(mat[i])):
                if i == j: # Zeroes the diagonal
                    adj[i][j] = 0
                elif mat[i][j] >= cutoff:
                    if map_type == HARD:
                        adj[i][j] = 1
                    else:
                        adj[i][j] = mat[i][j]
                else:
                    adj[i][j] = 0
    elif map_type == SOFT:
        for i in range(0,len(mat)):
            for j in range(0,len(mat[i])):
                if i == j: # Zeroes diagonal
                    adj[i][j] = 0
                else:
                    adj[i][j] = np.power((mat[i][j]+1)*0.5,beta)
    else: # This shouldn't happen
        raise TypeError("You chose an invalid map_type entry")
    if map_type == HARD:
        return np.array(adj).astype(int)
    return np.array(adj)

def network_measures(mat,weighted=False,limited=False,gamma=1.0):
    """
    Returns a dictionary of network statistics, including but not limited to: node
    degree, clustering coefficient, assortativity, local/global efficiency,
    modularity. Also returns centrality measures such as PageRank, betweenness
    centrality, etc.

    Input:
        mat: n x n 2D adjacency matrix, preferably with 0's on diagonal. Range [0,1] (non-negative weights).
        weighted: Default False; modifies whichever calculations to indicate that the matrix has weighted vertices
        gamma (opt): passed to modularity/community structure methods, lower than 1 preferentially searches for larger modules
        limited: only run/return a limited subset of measures to speed up operation, particularly as some of the excluded
            measures may not be well-defined for our graphs

    Output:
        dictionary with the statistic names as keys, results/arrays as values
    """
    debug_timing = False
    if debug_timing:
        import time
        currtime = time.clock()
    measures = {}
    measures["num_edges"] = np.count_nonzero(mat)
    measures["num_vertices"] = len(mat) # Lazy, but even a vertex with no connections should be counted
    # degree
    measures["degree"] = bct.bct.degrees_und(mat)
    measures["mean_degree"] = np.mean(measures["degree"])
    if weighted:
        measures["strength"] = bct.bct.strengths_und(mat)
        measures["mean_strength"] = np.mean(measures["strength"])
    if debug_timing:
        newtime = time.clock()
        delta = newtime - currtime
        currtime = newtime
        print "degree",delta
    # clustering coeff
    measures["clustering_binary"] = bct.bct.clustering_coef_bu(mat)
    measures["mean_clustering_binary"] = np.mean(measures["clustering_binary"])
    if weighted:
        measures["clustering_weighted"] = bct.bct.clustering_coef_wu(mat)
        measures["mean_clustering_weighted"] = np.mean(measures["clustering_weighted"])
    if debug_timing:
        newtime = time.clock()
        delta = newtime - currtime
        currtime = newtime
        print "clustering coeff",delta
    # assortativity
    if not limited:
        measures["assortativity_binary"] = bct.bct.assortativity_bin(mat, 0)
        measures["mean_assortativity_binary"] = np.mean(measures["assortativity_binary"])
        if weighted:
            measures["assortativity_weighted"] = bct.bct.assortativity_wei(mat, 0)
            measures["mean_assortativity_weighted"] = np.mean(measures["assortativity_weighted"])
        if debug_timing:
            newtime = time.clock()
            delta = newtime - currtime
            currtime = newtime
            print "assortativity",delta
    # characteristic path length
    # global efficiency
    if not limited:
        if weighted:
            distance = bct.bct.distance_bin(mat) # Of note, weighted distance not currently used because algorithm requires length matrix (smaller numbers indicate stronger connection) as opposed to weights
        else:
            distance = bct.bct.distance_bin(mat)
        cp = bct.bct.charpath(distance) # returns charpath,efficiency,ecc,radius,diameter
        measures["charpath"]=cp[0]
        measures["global_efficiency"] = cp[1]
        measures["eccentricity"]=cp[2]
        measures["radius"]=cp[3]
        measures["diameter"]=cp[4]
        if debug_timing:
            newtime = time.clock()
            delta = newtime - currtime
            currtime = newtime
            print "characteristic path length",delta
    # Local efficiency - inverse of char path length; cannot use weighted measure as it requests a weighted distance matrix
    measures["local_efficiency"]=bct.bct.efficiency_bin(mat,local=True)
    measures["mean_local_efficiency"] = np.mean(measures["local_efficiency"])
    if debug_timing:
        newtime = time.clock()
        delta = newtime - currtime
        currtime = newtime
        print "Local efficiency ",delta
    # modularity
    if not limited:
        mod = bct.bct.modularity_louvain_und(mat,gamma)
        measures["modularity"] = mod[1]
        measures["community_structure"] = mod[0]
        if debug_timing:
            newtime = time.clock()
            delta = newtime - currtime
            currtime = newtime
            print "modularity",delta
    # Giant component
    if not limited:
        measures["giant_component"] = np.max(bct.bct.get_components(mat)[1])
        if debug_timing:
            newtime = time.clock()
            delta = newtime - currtime
            currtime = newtime
            print "Giant component",delta
    # TODO: Ratio of mean clustering coefficient to mean clustering coefficient in randomly wired network with same degree distribution (C/C_{ran})
    # Ratio of C/C_ran to L/L_ran (L_ran is characteristic path length of randomly wired network with same degree distribution)

    # betweenness
    measures["betweenness_binary"] = bct.bct.betweenness_bin(mat)
    #if weighted: # NOT USED currently because the function needs a "connection-length" matrix, which would be inverted from our current weighted matrix
        #measures["betweenness_weighted"] = bct.bct.betweenness_weighted(mat)
    if debug_timing:
        newtime = time.clock()
        delta = newtime - currtime
        currtime = newtime
        print "betweenness",delta
    # eigenvector
    measures["eigenvector_centrality_und"] = bct.bct.eigenvector_centrality_und(mat)
    if debug_timing:
        newtime = time.clock()
        delta = newtime - currtime
        currtime = newtime
        print "eigenvector",delta
    # Pagerank
    measures["pagerank"] = bct.bct.pagerank_centrality(mat,0.85) # MAGIC NUMBER WARNING; default pagerank dampening factor
    if debug_timing:
        newtime = time.clock()
        delta = newtime - currtime
        currtime = newtime
        print "Pagerank",delta
    return measures

def compare_networks(arr_of_measure_dicts,names,list_of_measures = None):
    """
    Crude method to print a table comparing the above-generated measures between
    graphs. Gives good basic idea. Does not print lists, only singular values
    """
    print "{:30s}".format("measure"),
    for i in names:
        print "{:15s}".format(i[:14]),
    print ""

    if list_of_measures == None:
        set_of_all_measures = set()
        for s in arr_of_measure_dicts:
            for k,v in s.iteritems():
                if k not in set_of_all_measures:
                    set_of_all_measures.add(k)
        list_of_measures = sorted(list(set_of_all_measures))

    for k in list_of_measures:
        print "{:30s}".format(k[:29]),
        for i in range(0,len(arr_of_measure_dicts)):
            if k in arr_of_measure_dicts[i]:
                v = arr_of_measure_dicts[i][k]
                if type(v) != np.ndarray:
                    print "{:15s}".format(str(v)[:14]),
                else:
                    print "{:15s}".format("len([])={:d}".format(len(v))),
            else:
                print "{:15s}".format(" "),
        print ""

def print_nodal_measures(measure_dict):
    arr_names = []
    arr = []
    size = measure_dict["num_vertices"]
    for k,v in measure_dict.iteritems():
        if type(v) == np.ndarray and len(v) == size:
            arr.append(v)
            arr_names.append(k)
    print "nodes,"+','.join(arr_names)
    for i in range(0,size):
        print str(i)+','+','.join([str(measure[i]) for measure in arr])


def get_ROI_list(loadfile):
    f = open(loadfile,'r')
    roi = [line.strip().split('\t') for line in f]
    f.close()
    return roi

def write_ROI_node_file(loadfile,writefile,coloring_list,poi_list):
    """
    Writes an ROI node file for BrainNetViewer using a template file, rewriting
    coloring based on the community structure given in coloring_list, and increases
    the size of the POI list
    Input:
        loadfile: name of file to input, contains N entries
        writefile: name of file to output (will overwrite prior files there!)
        coloring_list: list of length N, each entry numbered according to group
        poi_list: list of names of ROIs which we will look up on ROI table to find appropriate ROI

        todo: Eventual goal is to pass optional poi_index directly to this, bypassing the named array
    """
    roi = get_ROI_list(loadfile)
    for i in range(0,len(roi)):
        roi[i][3] = int(coloring_list[i])
    poi_index = [mni.find_closest(mni.common_roi[i],roi)[0] for i in poi_list]
    for p in poi_index:
        roi[p][4]=4 # Default size is 2
    f = open(writefile,'w')
    for i in range(0,len(roi)):
        for item in roi[i]:
            f.write(str(item)+'\t')
        f.write('\n')
    f.close()

def network_measures_helper_generator(args):
    """
        Used for multiprocessing; optional dictionary with keys 'weighted', 'limited', 'gamma'
        Returns:
            partial function with above args filled in, only requiring mat
    """
    return partial(network_measures,**args)

def parallel_function(f):
    def easy_parallelize(f, sequence,pool_size=8):
        """ assumes f takes sequence as input, easy w/ Python's scope """
        from multiprocessing import Pool
        pool = Pool(processes=pool_size) # depends on available cores
        result = pool.map(f, sequence) # for i in sequence: result[i] = f(i)
        cleaned = [x for x in result if not x is None] # getting results
        cleaned = np.asarray(cleaned)
        pool.close() # not optimal! but easy
        pool.join()
        return cleaned
    from functools import partial
    return partial(easy_parallelize, f)


if __name__ == '__main__':
    mat_control,mask_control = load_patients([control_folder+'/'+i+'/'+filename for i in control]) # will be mat[x][y][i] where i is trial, x and y represent the correlation adjacency matrix for that patient on the Power 264 node setup
    adj_control = average_patients(mat_control,mask_control)
    G_control,sparsed_adj_control = create_graph(adj_control,95,weighted=False)

    mat_population,mask_population = load_patients([population_folder+'/'+i+'/'+filename for i in population]) # will be mat[x][y][i] where i is trial, x and y represent the correlation adjacency matrix for that patient on the Power 264 node setup
    adj_population = average_patients(mat_population,mask_population) # adj_population[x][y] represents averaged adjacency matrix for all trials
    G_population,sparsed_adj_population = create_graph(adj_population,95,weighted=False)

    #ig.summary(G_control)

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
    clusters_population = G_population.community_infomap()
    membership_population = clusters_population.membership
    # ig.plot(cluster, vertex_label=range(0,len(adj_control)),vertex_label_size=8,bbox=[1000,1000]) # PLOT community clusters

    #Community comparators 
    #f1 = open('power_communities.txt','r') # No longer can use because of new ROI list
    #power = [int(f.strip()) for f in f1.read().split('\n') if len(f)>0]
    #f1.close()
    #print "NMI score (control to Power): ",ig.clustering.compare_communities(power,membership_control,method="nmi") # nmi, vi, etc
    print "NMI score (population to control): ",ig.clustering.compare_communities(membership_population,membership_control,method="nmi") # nmi, vi, etc

    f = open('ROI_nodes_new_v2.node','r')
    roi = [line.strip().split('\t') for line in f]
    f.close()
    for i in range(0,len(membership_control)):
        roi[i][3] = membership_control[i]
    poi_index,dist = mni.find_closest(mni.common_roi['l_amygdala'],roi)

    
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
    #print_nodal_measures(network_measures(adj_control, weighted = True))

    colored_mat_control = [[membership_control[j] if membership_control[j]== membership_control[i] else -2 for j in range(0,len(adj_control[i])) ] for i in range(0,len(adj_control))]
    colored_mat_population = [[membership_population[j] if membership_population[j] == membership_population[i] else -2 for j in range(0,len(adj_population[i])) ] for i in range(0,len(adj_population))]
    pylab.subplot(2,2,1)
    draw_corr_matrix(adj_control,show=False)
    pylab.subplot(2,2,2)
    draw_corr_matrix(colored_mat_control,show=False)
    pylab.subplot(2,2,3)
    draw_corr_matrix(sort2d(adj_control,membership_control),show=False)
    pylab.subplot(2,2,4)
    draw_corr_matrix(sort2d(colored_mat_control,membership_control),show=False)
    pylab.show()
    #write_ROI_node_file('ROI_nodes_new_v2.node','my_ROI_new_v2.node',membership_control,['l_amygdala'])
    
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

