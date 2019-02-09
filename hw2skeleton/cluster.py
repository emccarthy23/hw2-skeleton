from .utils import Atom, Residue, ActiveSite
from .io import generate_distance_types
def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    similarity = 0.0
    list_size = 0
    global_counter = 0
    distance_types = generate_distance_types()
    """
	My similarity metric is based off of the algorithm PocketMatch (Yeturu. 2008. BMC Bioinformatics). 
	It uses the distance between three representative points (the alpha carbon, the beta carbon, and 
	the centroid of the side chain atoms) in each residue to create a set of distances to classify each active site. 
	The distances are sorted into 90 groups based on the atom type and group type for the two residues within 
	biologically relevant amino acids groups (e.g. Group-0:Hydrophobic; Group-1:Positive charge; Group-2:Negative charge;
	Group-3:Aromatic; Group-4:Polar). 
	The similarity score for two active sites is calculated by calculating the net average of how many distances 
	within distance type for two residues are within a threshold of 0.5 A (this distance cutoff was recommended 
	in the original algorithm based on the average amount of error in a PDB structure). 
	For my similarity metric, I converted this similarity score to a dissimilarity score (1 – similarity score) 
	so the values would be interpretable as a “distance”.
	Step 1: Assign residues in each active site to the following groups (This is coded into read_active_site). Group-0:Hydrophobic; Group-1:Positive charge; Group-2:Negative charge; Group-3:Aromatic; Group-4:Polar
	
	Step 2: Compute the Centroid for each residue in each active site (This is coded into the read_active_site). The centroid is the centroid coordinate for the side chain atoms (not including CB)
	
	Step 3: Compute distances between all residues within one active site (This is coded into the read_active_site). Compute distance between CA, CB, and Centroid for each residue in each active site. Bin distances in ascending order by "type" with each groupA,groupB,atomA,atomB being a type (so 90 groups since there are 5 groups and 3 atom types)
	
	Step 4: Sort each distance list in ascending order (This is coded into the read_active_site).
	
	Step 5: Align distances for two active sites. For site_a and site_b, determine for each type of distance how many distances are within a distance of 0.5 Angstroms.
	Step 6: Score alignment. Create score for alignment by taking sum of all distances that aligned divided by sum of the max of the length for each distance list between the two residues. Score of 1 is perfect alignment and score of 0 is no alignment
	Step 7: Change similarity score into dissimilarity score. Convert score to 1- score so that a score of 0 is for the same residues and a score of 1 is for residues with zero aligning distances

    """
    
    tau = 0.5
    for key in distance_types:
        #Subset distance lists to compare
        S1 = site_a.distances[key]
        S2 = site_b.distances[key]
        #Initialize parameters
        i = 0
        m = len(S1)
        j = 0
        n = len(S2)
        counter = 0
        #Initialize parameters
        if len(S1) != 0 and len(S2) != 0:
            while (i < m) and (j < n):
                if abs(S1[i] - S2[j]) <= tau:
                    i = i + 1
                    j = j+1
                    counter = counter + 1
                else:
                    if S1[i] < S2[j]:
                        i = i + 1
                    else:
                        j = j + 1
            global_counter = global_counter + counter
        list_size = list_size + max(m,n)
    if list_size != 0:
        similarity = 1 - global_counter/list_size

    return similarity

import random
def cluster_by_partitioning(active_sites,number_clusters):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Number of clusters
    k = number_clusters
    cluster_dict = {}
    # Initialize cluster labels
    # Choose random centroids
    C = random.sample(range(len(active_sites)), k)
    for i in range(k):
        cluster_dict[i] = [C[i]]

    done = False
    # Loop will run till the error becomes zero
    while not done:
        # Assigning each value to the cluster with the closest centroid
        for site_id in (set(list(range(len(active_sites)))) - set(C)):
            distances = []
            for centroid in C:
                distances.append(d[site_id,centroid])
            cluster = np.argmin(distances)
            cluster_dict[cluster] = cluster_dict[cluster] + [site_id]
        # Storing the old centroid values
        C_old = C
        # For each cluster, find the sum of the distances to each element of the cluster to all the other elements.
        #Take the element with the minimum sum as the new centroid
        for cluster_id in range(k):
            cluster_indices = cluster_dict[cluster_id]
            #if cluster_indices == []:
            #    return [], []
            #else:
            new_centroid = cluster_indices[0]
            dist_sum = len(active_sites) -1
            for index_1 in cluster_indices:
                index_dist = []
                for index_2 in cluster_indices:
                    index_dist.append(d[index_1,index_2])
                if np.sum(index_dist) < dist_sum:
                    dist_sum = np.sum(index_dist)
                    new_centroid = index_1
            C[cluster_id] = new_centroid
            #print(C)
        #End if the centroids didn't update
        if C == C_old:
            done = True
    #print(C)
    return list(cluster_dict.values())

def create_dendogram(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    #Compute a distance matrix for all the sites
    n = len(active_sites)
    d = np.zeros((n,n))
             
    for i in range(n):
        for j in range(i,n):
            d[i,j] = compute_similarity(active_sites_list[i], active_sites_list[j])
            d[j,i] = d[i,j]
   
    #Converts distance matrix to an upper triangular matrix
    y = pdist(d)
    """
    Computes a dendogram for the sites using single linkage as the cluster combining metric
    This is an agglomerative clustering meaning it is using a bottom up approach
    First every site is assigned to its own cluster. Then each cluster is combined with the cluster it is cloesest to
    based on the single linkage distance metric. For thi smetric, the distance between two clusters is defined as
    the shortest distance between an element in the first cluster and an element in the second cluster. This process
    is iteratively repested until there is one cluster
    """
    Z = single(y)
    return Z
    
def cluster_hierarchically(active_sites,distance):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    #Compute a distance matrix for all the sites
    Z = create_dendogram(active_sites)
    n = len(active_sites)
    #This flattens the denogram that was created in the previous step. So clusters are created from the dendogram with the
    #criterion that "t" is the max distance in a cluster 
    cluster_assignments = fcluster(Z, t=distance, criterion='distance')
    #Store the sites in a list of lists where each site is in a list with the other sites in the same cluster
    number_of_clusters = len(np.unique(cluster_assignments))
    clusters = [[]]*number_of_clusters
    for i in range(number_of_clusters):
        clusters[i] = np.array(range(n))[cluster_assignments == (i+1)]
    return clusters



