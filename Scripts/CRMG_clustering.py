import numpy as np
import sys, os
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import glob
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from sklearn.cluster import DBSCAN
# split KH, RRM
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy import sparse

# Defines a human readable name for a collection of pre-defined sublclades
# Sub clades summarize sets of species
def bestmatchedcladename(klist):
    clades = [['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm', 'Algae', 'Plant', 'Flowering_plant', 'Protists', 'Plasmodium', 'Plant_pathogen', 'Heterokonts'], ['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm', 'Algae', 'Plant', 'Flowering_plant', 'Protists'], ['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm', 'Algae', 'Plant', 'Flowering_plant'], ['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm'], ['Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm'], ['Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate'], [ 'Vertebrate', 'Bird', 'Mammal', 'Primate'], ['Mammal', 'Primate'], ['Insect', 'Fly', 'Worm'], ['Insect', 'Fly'], ['Algae', 'Plant', 'Flowering_plant'], ['Plant', 'Flowering_plant'], [ 'Plant_pathogen', 'Heterokonts'], ['Plasmodium', 'Plant_pathogen', 'Heterokonts']]
    cladenames = ['LECA', 'LECA-some', 'Fungi+Metazoan+Plant', 'Fungi+Metazoan','Metazoan','Vertebrates','Land-vertebrates','Mammals', 'Protostomia', 'Insects', 'Green plants', 'Land plants', 'Stramenopiles', 'Alveolata']
    
    if len(klist)>1:
        fraction = []
        for cla in clades:
            frac = float(len(np.intersect1d(cla,klist)))/float(len(np.union1d(cla, klist)))
            fraction.append(frac)
        return cladenames[np.argmax(fraction)]
    else:
        return klist[0]

def assignclade(cluster):
    assignedclade = []
    for c, clust in enumerate(cluster):
        if np.sum(clust>0) > 1:
            assclust = evodistance[clust>0][0, :]<=np.amax(evodistance[clust>0][:, clust>0])
            asskings = np.unique(kings[assclust])
            assignedclade.append(bestmatchedcladename(asskings))
        else:
            assignedclade.append(species[clust>0][0])
    return assignedclade


# control function that finds string in array of strings, string might be a substring in the strings in array
def findstring(array, string, case = False):
    if case:
        string = string.upper()
    for a in array:
        if case:
            a = a.upper()
        if string in a:
            return True
    return False


                       

# compute homolog matrix, 1-1 ortholog matrix and triangle 1-1 ortholog matrix pairs
# these matrices contain additional information about the sequence identity relationship between species
def homologmatrix(idmat, pspecies, testortholog = False, testtriangle = False, maskmat = None, include_self = False, min_evidence = 30):
    uspec = np.unique(pspecies)
    hmatrix = np.zeros(np.shape(idmat), dtype = np.int8)
    minid = np.amin(idmat)
    # assign best neighbor for each protein to a protein in all other species
    for ui in uspec:
        uiloc = np.where(pspecies == ui)[0]
        for uj in uspec:
            if ui == uj:
                hmatrix[uiloc, uiloc] = 1
            else:
                ujloc = np.where(pspecies == uj)[0]
                midmat = idmat[uiloc][:, ujloc]
                if maskmat is not None:
                    midmat[~maskmat[uiloc][:,ujloc]] = minid
                dneigh = np.argmax(midmat ,axis = 1)
                hmatrix[uiloc, ujloc[dneigh]] = 1
    #defines minimum sequence identity between two RBPs to be homologs
    if min_evidence is not None:
        print('Homolog connections before 30% cutoff', np.sum(hmatrix))
        hmatrix = hmatrix * (idmat > min_evidence ).astype(np.int8)
        print('after', np.sum(hmatrix))
    # if proteins point at each other, they can be orthologs
    if testortholog:
        omatrix = hmatrix * hmatrix.T == 1
        
    else:
        omatrix = None
    # if proteins are orthologs and also share another ortholog, they form a triangle, which is a strong indicator of common ancestor
    if testtriangle:
        A = sparse.csr_matrix(omatrix.astype(int))
        tmatrix = A.dot(A)
        tmatrix = tmatrix.toarray()*omatrix.astype(int)
        tmatrix = tmatrix >= 3
        np.fill_diagonal(tmatrix, True)
        
    else:
        tmatrix = None
    # whenever one protein is closest to another, they could be homologs
    if include_self:
        for ui in uspec:
            uiloc = np.where(pspecies == ui)[0]
            midmat = np.argsort(idmat[uiloc][:, uiloc], axis = 1)[:,-2]
            hmatrix[uiloc, uiloc[midmat]] = 1
            hmatrix[uiloc[midmat], uiloc] = 1
    hmatrix = (hmatrix + hmatrix.T) >= 1
    
    return hmatrix, omatrix, tmatrix


def tree_guided_clustering(treedist, leaftype, valuemat, cutoff, cmask = None):
    # Inputs: 
        # treedist: matrix that contains the age distances between pairs of species
        # leaftype: species of proteins given to the clustering
        # valuemat: JPLE distances of proteins in the clustering
        # cutoff: JPLE distance cutoff for clusters
    # unique species in the set of rbps to cluster
    unspecies = np.unique(leaftype)
    # Determine all species that are in the clade, formed by the set of species. F.e. if one species was missed by metazoan proteins, include this species
    included_species = np.where(treedist[unspecies[0]]<= np.amax(treedist[unspecies][:,unspecies]))[0]
    # from this list of species in the clade, determine all the brachpoints towards the most common ancestor
    evo_steps = np.unique(treedist[included_species][:,included_species])[1:] # look at all evolutionary distances sequentially, excluding 0 
    # (IMPORTANT: distances between species of different clades in treedist need to be non-redundant, otherwise distance branches are combined)
    
    # assign intial cluster numbers: at the beginning, every RBP is its own cluster
    pot_clusters = np.arange(len(leaftype)) # this array one will be reduced with the valuemat to reflect the number of clusters
    # this array counts how many times no merge was performed when moving up the tree. If no merge occured after one two branches, this cluster cannot be merged anymore.
    descendant = np.ones(np.size(pot_clusters))
    
    iscluster = np.copy(pot_clusters) # this array will stay the same size to assign clusters to individual proteins
    # To average the JPLE distance of two branches over the entities of the two clusters we have to keep track how many rpbs are in the cluster
    counts = np.ones(np.shape(valuemat))
    # the similarity between themselves is not accounted for and therefore the counts between are set to 0
    np.fill_diagonal(counts, 0)
    # start with the most recent branchpoint and move to the most distant.
    for l, lt in enumerate(evo_steps):
        # select all species that all under this brachpoint
        lt = np.unique(included_species[np.where(treedist[included_species][:,included_species] == lt)[0]])
        # select all proteins from these species
        potp = np.where(np.isin(leaftype, lt))[0]
        
        '''
        # control segment
        if cmask is not None:
            oneis = 0
            if 'P086417_1.95d' in latentnames[cmask[potp]]:
                oneis += 1
                piii = iscluster[latentnames[cmask] == 'P086417_1.95d'][0]
                pjjj = iscluster[latentnames[cmask] == 'P103080_1.95d'][0]
                print(realprotnames[cmask][iscluster == piii], realprotnames[cmask][iscluster == pjjj])
                print(valuemat[pot_clusters == piii, pot_clusters == pjjj]/counts[pot_clusters == piii, pot_clusters == pjjj])
                
                print(l, 'P086417_1.95d', descendant[np.isin(pot_clusters,piii)], descendant[np.isin(pot_clusters,pjjj)])
                
            if 'P103080_1.95d' in latentnames[cmask[potp]]:
                oneis += 1
                piii = iscluster[latentnames[cmask] == 'P086417_1.95d'][0]
                pjjj = iscluster[latentnames[cmask] == 'P103080_1.95d'][0]
                print(realprotnames[cmask][iscluster == pjjj], realprotnames[cmask][iscluster == piii])
                print(valuemat[pot_clusters == pjjj, pot_clusters == piii]/counts[pot_clusters == pjjj, pot_clusters == piii])
                print(l, 'P103080_1.95d', descendant[np.isin(pot_clusters,pjjj)], descendant[np.isin(pot_clusters,piii)])                
        '''
        
        # increase the number of descendant by 1 for all the clusters that these proteins are representing.
        descendant[np.isin(pot_clusters, iscluster[potp])] += 1
        # Perform iterative process to avoid ambigouties through multiple clusters that can be merged at once. Always merge closest clusters first and then recompute the merging statistics
        while True:
            # Get positions of clusters in that species in the count and valuemat
            potc = np.where(np.isin(pot_clusters,iscluster[potp]))[0]
            # Only include clusters that fullfil direct decendant rule from current brachpoint, i.e where descendant has not been increased twice
            potc = potc[descendant[potc] <= 3]
            # stop iteration if there is no more two cluster to potentially merge
            if len(potc) <=1:
                break
            # get the valuemat for these clusters
            dvmat = valuemat[potc][:,potc]
            # transform value mat, which is a sum of all values between all elements in that matrix into a per pair mean matrix
            dvmat = dvmat/counts[potc][:,potc]
            # set the diagonal of the valuemat to the max so that we can use np.argmin to find the most similar pairs
            np.fill_diagonal(dvmat, cutoff +1)
            # set everything to max that's not within the JPLE cutoff
            dvmat[dvmat >= cutoff] = cutoff +1
            if (dvmat <= cutoff).any(): # check if any cluster pair can be merged
                # find closest pair to be united, closepair contains the coordinates
                closepair = [int(np.argmin(dvmat)/len(dvmat)), np.argmin(dvmat)%len(dvmat)]
                # transform into coordinates in the orignal value matrix and arrays
                potc = potc[closepair]
                # unite pair in counts, valuemat and pot_clusters
                iscluster[iscluster == pot_clusters[potc[0]]] = pot_clusters[potc[1]]
                pot_clusters[potc[0]] = pot_clusters[potc[1]]
                # if cluster merged, set descendant to 1 so that cluster can be considered for merging in three more steps
                # if all the cluster were actually from the same species, then set them to 2 because then they actually missed this branchpoint and only get one more change to be merged.
                if len(np.unique(leaftype[np.isin(iscluster, pot_clusters[potc])])) == 1:
                    descendant[potc[0]] = 2
                else:
                    descendant[potc[0]] = 1
                # delete the merged cluster
                pot_clusters = np.delete(pot_clusters, potc[1])
                descendant = np.delete(descendant, potc[1])
                # add counts from row and columns to each other and remove rows and colums from merged cluster. Adding up both dimensions will result in the number of pairs that are compared between each other. So that the sum divided by the number will be the average JPLE distance between pairs in that group.
                counts[potc[0]] += counts[potc[1]]
                counts = np.delete(counts, potc[1], axis = 0)
                counts[:, potc[0]] += counts[:, potc[1]]
                counts = np.delete(counts, potc[1], axis = 1)
                
                valuemat[potc[0]] += valuemat[potc[1]]
                valuemat = np.delete(valuemat, potc[1], axis = 0)
                valuemat[:, potc[0]] += valuemat[:, potc[1]]
                valuemat = np.delete(valuemat, potc[1], axis = 1)
                
            else:
                break
    # sort clusters by size
    ucluster, uclustersize = np.unique(iscluster, return_counts = True)
    ucluster = ucluster[-np.argsort(uclustersize)]
    clusterout = np.copy(iscluster)
    # rename the clusters with new integers and assign the mean jple between objects to them as a confidence
    clusterval = np.zeros(len(ucluster))
    for c, cl in enumerate(ucluster):
        clusterval[c] = valuemat[list(pot_clusters).index(cl), list(pot_clusters).index(cl)]/max(1,counts[list(pot_clusters).index(cl),list(pot_clusters).index(cl)])
        clusterout[iscluster == cl] = c
    cluster = clusterout
    return cluster, clusterval


# clustering which combines all nodes connected by a network matrix with True and False, True indicating that there is a connection
def clusterall(fmat):
    cluster = -np.ones(len(fmat))
    c = 0
    cint = [0]
    while True:
        expand = np.where(np.sum(fmat[cint], axis = 0)>0)[0]
        if len(expand) > len(cint):
            cint = np.copy(expand)
        else:
            cluster[cint] = c
            c+=1
            cint = np.where(cluster == -1)[0]
            if len(cint) == 0:
                break
            else:
                cint = [np.where(cluster == -1)[0][0]]
    return cluster


# function uses idmat to cluster with scipy linkage, using method as linkage
def conserved_feature_group_clustering(hidmat, flatdist, protspecies, confcut, maskmat = None):
    # hidmat: boolean matrix that has true between entries that are connected by sequence homology
    # flatdist: cosine distance from JPLE
    # protspecies: species id for each protein
    # confcut: JPLE distance cutoff for clusters
    # maskmat: additional boolean contraint matrix that accounts for domain composition for example.
    
    if maskmat is None:
        maskmat = hidmat
    
    # combine all RBRs that are homolog connected and have similar binding according to flatdist
    cluster = clusterall(hidmat*(flatdist<=confcut) * maskmat)
    cls, cln = np.unique(cluster, return_counts = True)
    print('Initial homolog clusters', len(np.unique(cluster)))
    ccsizes, ccnsizes = np.unique(cln, return_counts = True)
    print('Size distribution', '\n', ccsizes, '\n', ccnsizes)
    
    # resort clusters by size
    ncluster = -np.ones(len(cluster))
    cls = cls[np.argsort(-cln)]
    for c, cl in enumerate(cls):
        ncluster[cluster == cl] = c
    cluster = ncluster
    maxcl = np.amax(cluster)

    for c in range(len(cls)):
        cl = cls[c]
        pstuff = False
        cmask = np.where(cluster == cl)[0]
        
        # Control segment to see if algorith works
        #isinthe = False
        #if 'P086417_1.95d' in latentnames[cmask]:
            #print(realprotnames[cmask], len(cmask))
            #isinthe = True
        
        if len(cmask) > 1:
            cmat = flatdist[cmask][:,cmask]
            # refine clusters based on their phylogenetic composition and average distance to each other
            nc, ncdist = tree_guided_clustering(evodistance, protspecies[cmask], cmat, confcut, cmask = cmask)
            for n in np.unique(nc)[1:]:
                maxcl +=1
                cluster[cmask[nc == n]] = maxcl
                
                ### second part of controls segment
                #if isinthe:
                    #if 'P103080_1.95d' in latentnames[cmask[nc == n]] or 'P086417_1.95d' in latentnames[cmask[nc == n]]:
                        #print(n, realprotnames[cmask[nc == n]])
                
    cls, cln = np.unique(cluster, return_counts = True)
    
    ncluster = -np.ones(len(cluster))
    cls = cls[np.argsort(-cln)]
    for c, cl in enumerate(cls):
        ncluster[cluster == cl] = c
    cluster = ncluster.astype(int)
    maxcl = np.amax(cluster)
    return cluster

# function splits linkage_matrix into clusters at a cutlink cutoff, or a alldistabove matrix that is True for every link that is above the cluster
def splittree(linkage_matrix, cutlink, alldistabove = None):
    le = len(linkage_matrix)+1
    if alldistabove is None:
        alldistabove = linkage_matrix[:,2]>=cutlink
    
    cluster = -np.ones(le, dtype = int)
    clusterdist = np.zeros(le)
    # takes quite long for single linkage, should be changed into top-down algorithm
    if np.sum(alldistabove) > 0:
        #for each leaf follow tree until hit cutoff and assign branch number to cluster
        for i in range(le):
            stage = np.copy(i)
            ite = 0
            while True:
                nstage = np.where(linkage_matrix[:, [0,1]] == stage)[0][0]
                if alldistabove[nstage]:
                    cluster[i] = stage
                    if ite == 0:
                        clusterdist[i] = 0
                    else:
                        clusterdist[i] = linkage_matrix[stage-le, 2]
                    break
                else:
                    stage = nstage + le
                    ite +=1
    # sort clusters
    ucluster, uclustersize = np.unique(cluster, return_counts = True)
    ucluster = ucluster[-np.argsort(uclustersize)]
    clusterout = np.copy(cluster)
    for c, cl in enumerate(ucluster):
        clusterout[cluster == cl] = c
    cluster = clusterout
    return cluster, clusterdist

# unite clusters within species if their assigned clusters share mean less than cutoff
# to correct for paralog clusters that weren't connected because of paralog 

# Algorithm: 
    # Find potential proteins within the same species that may have been forgotten in the origianl clustering because their homolog network does not interact with each other since they are paralogs.
    # Find proteins that are part of a cluster that is represented in the same or a subset of the species of the larger cluster and check if the two clusters can be combined
    # Check their mergeability through agglomerative clustering
        # If clusters or parts of both clusters are merged in agglomerative clustering, they can be merged.
def combine_paralog_groups(cluster, flatdist, protspecies, maskmat):
    # cluster contains currect cluster assignment of proteins
    # maskmat establishes which proteins can be combined based on their domain structure
    cls, cln = np.unique(cluster, return_counts = True)
    for c, cl in enumerate(cls):
        smask = np.where(cluster == cl)[0]
        # species that are in the current cluster
        potspecies = np.unique(protspecies[smask])
        if len(potspecies) > 0:
            potprot = []
            for ps in potspecies: # iterate over all species in the cluster
                spmask = smask[protspecies[smask] == ps]
                canhavelink = flatdist[spmask] # distmatrix for the members of the cluster that are from that species to all other proteins
                canhavelink[~maskmat[spmask]] = confcut+ 0.1 # set the distances of proteins with false domain structure higher than cutoff
                havelink = np.where((np.sum(canhavelink <= confcut,axis = 0)>0)&(protspecies == ps))[0] # determine all proteins from the same species that have a link to any of the proteins in that group
                potprot.append(havelink)
            potclust = cluster[np.concatenate(potprot)]
            potclust = potclust[potclust!= cl] # determine the clusters from which these potential proteins come
            potclust = np.unique(potclust)
            if len(potclust) > 0:
                # only consider clusters whose species are covered 100% by potspecies
                keep = [] 
                for pcl in potclust:
                    pclspecies = np.unique(protspecies[cluster == pcl])
                    if np.isin(pclspecies, potspecies).all():
                        keep.append(pcl)
                potclust = np.array(keep)
            if len(potclust) > 0:
                potprot = np.where(np.isin(cluster, potclust))[0] # find all proteins in these potential clusters
                cmask = np.union1d(smask, potprot) # determine union between proteins from cluster and proteins from other clusters
                cmat = flatdist[cmask][:,cmask] # distance matrix with original proteins in cluster and potential paralogs
                linkage_matrix = linkage(cmat[np.triu_indices(len(cmat),1)], method='average')
                nc, ncdist = splittree(linkage_matrix, confcut) # make new clusters from extended distance matrix
                
                ncwith = nc[cluster[cmask] == cl] # new assinged clusters of the original cluster proteins
                nctotest, nctestcount = np.unique(cluster[cmask[np.isin(nc,ncwith)]], return_counts = True) # original clusters of all the proteins in the clusters that are now assigned to the original cluster proteins
                nctotest = nctotest[-np.argsort(nctestcount)] # sort these cluster by their size
                nctotest = nctotest[nctotest != cl] # only choose the ones that are not the original one. 
                if len(nctotest)>0:
                    cluster[np.isin(cluster,nctotest)] = cl # and assign them the to the original one.
                    
    cls, cln = np.unique(cluster, return_counts = True)
    ncluster = -np.ones(len(cluster))
    cls = cls[np.argsort(-cln)]
    for c, cl in enumerate(cls):
        ncluster[cluster == cl] = c
    cluster = ncluster.astype(int)
    maxcl = np.amax(cluster)
    return cluster

# assign the cluster an age and generate a matrix that contains the number of proteins per species in the cluster
def assign_clusterage(clusterid, protspecies, evodistance):
    cluster = []
    clusterage = []
    clusterassign = np.unique(clusterid)
    for c, cl in enumerate(clusterassign):
        clspecies = protspecies[clusterid == cl] # find all species in cluster
        clspecies, clspecnum = np.unique(clspecies, return_counts = True)
        clpat = np.zeros(len(species), dtype = np.int8)
        clpat[clspecies] = clspecnum
        cluster.append(clpat)
        clusterage.append(np.amax(evodistance[clspecies][:, clspecies]))
    cluster = np.array(cluster, dtype = np.int8)
    clusterage = np.array(clusterage)
    return cluster, clusterage


# False positives can be spotted as sequences that are not one-to-one orthologs to a member of the cluster and would result in younger common ancestor if removed
def cleanfpcluster(evidmat, evodistance, clusterage, cluster, clusterid, pspecies):
    cmax = np.amax(clusterid)
    c = 0
    newcluster = []

    while True:
        cl = np.copy(cluster[c])
        cmask = np.where(clusterid == c)[0]
        cspec = cl > 0
        coverage = np.sum(evodistance[cspec] <= clusterage[c], axis = 0) > 0
        repeat = True
        if np.sum(coverage) > 3:
            notonetoone = np.sum(evidmat[cmask][:,cmask], axis = 1) == 1
            espec = pspecies[cmask]
            uspec, usnum = np.unique(espec, return_counts = True)
            uspec = uspec[usnum == 1]
            #if any rbp has no 1-to-1 ortholog test if its removal changes common ancestor of cluster
            if notonetoone.any():
                arange = np.arange(len(species))
                for nto in np.where(notonetoone)[0]:
                    if espec[nto] in uspec:
                        dspec = cspec * ~np.isin(arange, espec[nto])
                        ncluage = np.amax(evodistance[dspec][:,dspec])
                        if ncluage < clusterage[c]:
                            cmax += 1
                            clusterid[cmask[nto]] = cmax
                            cluster = np.append(cluster,[np.isin(arange,espec[nto]).astype(int)] ,axis = 0)
                            clusterage = np.append(clusterage, [0])
                            cluster[c, espec[nto]] = 0
                            cspec = cluster[c] > 0
                            clusterage[c] = ncluage
                            repeat = False
                            break
        if repeat:
            c += 1
            
        if c == len(cluster) - 1:
            break
    return cluster, clusterage, clusterid



# false negatives are defined as clusters that share more than uniteval identity to another cluster
def cleanfncluster(idmat, clusterid, uniteval, cluster, fnidcut = 70):
    uniquecluster, uniqueclen = np.unique(clusterid, return_counts = True)
    uniquecluster, uniqueclen = uniquecluster[np.argsort(uniqueclen)], np.sort(uniqueclen)
    
    idmat = np.copy(idmat)
    cidmat = np.zeros((len(uniquecluster), len(uniquecluster)))
    np.fill_diagonal(idmat, 0)
    print('Compute cluster identity matrix ...')
    for c, cl in enumerate(uniquecluster):
        for d in range(c+1,len(uniquecluster)):
            dl = uniquecluster[d]
            cidmat[c,d] = cidmat[d,c] = np.mean(idmat[clusterid == cl][:, clusterid == dl])
    print('Merge clusters with ID larger than', fnidcut)
    while True:
        if not (cidmat >= fnidcut).any():
            break
        tounite = np.where(cidmat >= fnidcut)
        toi, toj = tounite[0][0], tounite[1][0]
        was, isnow = uniquecluster[toi], uniquecluster[toj]
        clusterid[clusterid == was] = isnow 
        cidmat[toj] = (uniqueclen[toi]*cidmat[toi] + uniqueclen[toj]*cidmat[toj])/(uniqueclen[toi]+uniqueclen[toj])
        cidmat[:,toj] = cidmat[toj]
        cidmat = np.delete(np.delete(cidmat, toi, axis = 0), toi, axis = 1)
        np.fill_diagonal(cidmat, 0)
        uniqueclen[toj] = uniqueclen[toi] + uniqueclen[toj]
        uniqueclen, uniquecluster = uniqueclen[uniquecluster != was], uniquecluster[uniquecluster != was]
        sortcls = np.argsort(uniqueclen)
        uniqueclen, uniquecluster, cidmat = uniqueclen[sortcls], uniquecluster[sortcls], cidmat[sortcls][:,sortcls]
    
    
    # reset cluster identifiers
    clusteridout = np.copy(clusterid)
    for c, cl in enumerate(uniquecluster[::-1]):
        clusteridout[clusterid == cl] = c
    clusterid = clusteridout
    
    return clusterid




# for each cluster, use sequence identity or homolog matrix to identify most likely common ancestor of cluster
# and determine their relationship
def determine_family(clusters, clusterage, evodistance, clusterid, simmat, trust = None, maskmat = None, maxin = 1, maxo = 0):
    
    familyage = np.zeros(len(clusters)) # age of sequence cluster from merging the two family members
    family = -np.ones(len(clusters)) # cluster id of family 
    famevidence = np.zeros(len(clusters)) # sequence identity between the family members
    family_relation = np.zeros(len(clusters), dtype = np.int8) # evolutionary relationship between the two family members
    
    print('Determine family relationships between conserved motif groups')
    for p in range(len(clusters)):
        if p % 250 == 0:
            print(p)
        testedprot = clusterid == p
        # mask all the species that have the specificity 
        currspec = np.where(cluster[p] > 0)[0]
        paralogage = clusterage[p]
        
        outsidespec = np.argsort(evodistance[currspec[0]])
        evosorted = evodistance[currspec[0]][outsidespec]
        
        # select all species within the clusters clade including the ones that the specificity has been lost
        insidespec = outsidespec[evosorted <= paralogage]

        # select all species outside the cluster
        outsidespec = outsidespec[evosorted > paralogage]
        evosorted = evosorted[evosorted > paralogage]
        # check if outside spec exists or if cluster fills all eukaryotes
        if len(outsidespec) > 0:
            # only look in species that are maximally two branchpoints away from cluster
            sortevo = np.unique(evosorted)
            maxout = min(maxo, len(sortevo)-1)
            outsidespec = outsidespec[evosorted<=sortevo[maxout]]
            
            potential_family = np.where(np.sum(clusters[:, outsidespec], axis = 1)>0)[0]
            # only consider clusters that are two brachpoints into the outside branch to avoid biases towards single RBP clusters
            fstages = []
            for mi in range(maxout+1):
                considerstages = evodistance[outsidespec]
                for co, cost in enumerate(considerstages):
                    if sortevo[mi] in cost:
                        cost[cost > sortevo[mi]] = sortevo[mi]
                        #print(mi, species[outsidespec[co]], sortevo[mi], np.unique(cost), np.unique(cost)[-1-maxin:])
                        fstages.append(np.unique(cost)[-1-maxin:])
            fstages = np.unique(np.concatenate(fstages))
            fstages = fstages[fstages != 0]
            fstages = np.unique(np.append(fstages, sortevo))
            potential_family = potential_family[np.isin(clusterage[potential_family], fstages)]
            
            ## add single species clusters if their cluster is only two away from brachpoint.
            potentialunique = outsidespec[(np.sort(evodistance[outsidespec], axis = 1)[:,1] == sortevo[0])]
            if maxin >=2:
                for mi in range(2,maxin+1):
                    potentialunique = np.unique(np.append(potentialunique, outsidespec[(np.sort(evodistance[outsidespec], axis = 1)[:,mi] == sortevo[0])]))
            if len(potentialunique) > 0:
                potential_family = np.append(potential_family, np.where((np.sum(cluster[:, potentialunique]>0, axis = 1)>=1) * (np.sum(cluster>0,axis=1)==1))[0])
            
            insidepotential = np.where(np.sum(clusters[:, insidespec], axis = 1)>0)[0]
            insidepotential = insidepotential[clusterage[insidepotential] >= sortevo[0]]
            potential_family = np.unique(np.append(potential_family, insidepotential))
            
            cores = ''.join(np.isin(np.arange(len(species)), insidespec).astype(int).astype(str))
            
            # compute identity between cluster and potential family 
            idtofam = np.zeros(len(potential_family))
            if maskmat is None:
                for f, pf in enumerate(potential_family):
                    pfprot = clusterid == pf
                    idtofam[f] = np.mean(simmat[testedprot][:,pfprot])
            else:
                for f, pf in enumerate(potential_family):
                    pfprot = clusterid == pf
                    idtofam[f] = np.mean(simmat[testedprot][:, pfprot][:, np.sum(maskmat[testedprot][:, pfprot],axis = 0)>0])
                    
            idtofam = np.array(idtofam)
            idtofam = np.nan_to_num(idtofam)
            relation = potential_family[np.argmax(idtofam)]
            relationid = idtofam[np.argmax(idtofam)]
            family[p] = relation
            famevidence[p] = relationid
            
            # compute species that are on branch of related cluster
            relationcore = evodistance[cluster[relation]>0][0] <= np.amax(evodistance[cluster[relation]>0][:, cluster[relation]>0])
            relationage = clusterage[relation]
            # check if cores overalap
            if len(np.intersect1d(np.where(relationcore)[0], insidespec)) == 0:
                # purple
                family_relation[p] = 3
                seqmask = np.amax(evodistance[np.append(np.where(relationcore)[0], insidespec)][:, np.append(np.where(relationcore)[0], insidespec)])
                familyage[p] = seqmask
            else:
                # red of brown
                # check if core of cluster overlaps with its realted cluster
                relationcluster = np.where(clusters[relation] > 0)[0]
                seqmask = np.amax(evodistance[np.append(np.where(relationcore)[0], insidespec)][:, np.append(np.where(relationcore)[0], insidespec)])
                familyage[p] = seqmask
                if len(np.intersect1d(relationcluster, insidespec)) == 0:
                    # red
                    family_relation[p] = 2
                else:
                    # brown
                    family_relation[p] = 1
        else:
            familyage[p] = np.amax(evodistance)
            famevidence[p] = 100.
    if trust is not None:
        notrust = famevidence < trust
        family_relation[notrust] = 4
        
    return family_relation, familyage, family, famevidence
        
        
def make_protdomain_arrays(protsimnames, maxdom = 5):
    # maximal number of domains that are considered to be in one protein sequence
    # generate two matrices that indicate the RBD class to which the RBPs belong
    # this one contains ones for x-1, x, and x+1 domain counts
    protdomains = np.zeros((2*maxdom,len(protsimnames)), dtype = np.int8)
    # contains only one for exact number of domains
    protdomains2 = np.zeros((2*maxdom,len(protsimnames)), dtype = np.int8)
    # generate readable names
    realprotnames = np.copy(protsimnames)
    # counts the number of domains
    domaincounts = np.zeros(len(protsimnames))
    
    # remove rbps with no info on domains
    remove = np.zeros(len(protsimnames)) == 0
    
    # generate arrays that determine which proteins can be compared to each other based on their domain composition
    # for example 1 rrm can be compared to 1 and 2 rrms but not 3
    for p, psname in enumerate(protsimnames):
        psname = psname.split('__')
        protsimnames[p] = psname[0]
        
        # Remove human RBPs that have no name
        if 'ENSG0' in psname[2] or ('KH' not in psname[-1] and 'RRM' not in psname[-1]): 
            remove[p] = False
        
        # check domain class is in name, then split domain names
        if '-' in psname[-1]:
            dtypes, dtnum = np.unique(psname[-1].split('-'), return_counts = True)
            dtypes = [dtypes[np.argmax(dtnum)]]
            dtnum = [np.amax(dtnum)]
        else:
            dtypes, dtnum = [psname[-1]], [1]
    
        # assign 1's to the entry which determines their potential interaction partners
        if 'RRM' in dtypes:
            lok = [min(max(0,x),maxdom-1) for x in range(-2+min(maxdom+1,dtnum[list(dtypes).index('RRM')]), 1+min(maxdom+1,dtnum[list(dtypes).index('RRM')]))]
            # arrays that define every domain type +-1 that can be similar to number of domains in the protein
            protdomains[lok,p] = 1
            #print(dtypes,dtnum, protdomains[:,p])
            # exact number of domains
            protdomains2[-1+min(maxdom,dtnum[list(dtypes).index('RRM')]), p] = 1
            #print(protdomains2[:,p])
            
        if 'KH' in dtypes:
            lok = [maxdom+min(max(0,x),maxdom-1) for x in range(-2+min(maxdom+1,dtnum[list(dtypes).index('KH')]), 1+min(maxdom+1,dtnum[list(dtypes).index('KH')]))]
            protdomains[lok,p] = 1
            protdomains2[maxdom-1+min(maxdom,dtnum[list(dtypes).index('KH')]), p] = 1
        
        domaincounts[p] = np.sum(dtnum)
        realprotnames[p] = psname[-2]+'__'+psname[-1]
 
    protsimnames = protsimnames[remove]
    realprotnames = realprotnames[remove]
    protdomains = protdomains[:, remove]
    protdomains2 = protdomains2[:, remove]
    domaincounts = domaincounts[remove]
    return protsimnames, realprotnames, domaincounts, protdomains, protdomains2, remove
     
        
def generate_cluster_stats(clusterassign, clusterid, realprotnames, protspecies, species):        
    speciesrank = ['Homo_sapiens', 'Mus_musculus', 'Macaca_mulatta', 'Gallus_gallus', 'Danio_rerio', 'Xenopus_tropicalis', 'Saccharomyces_cerevisiae', 'Schizosaccharomyces_pombe', 'Drosophila_melanogaster', 'Drosophila_ananassae', 'Caenorhabditis_elegans', 'Caenorhabditis_briggsae', 'Arabidopsis_thaliana', 'Arabidopsis_lyrata', 'Zea_mays', 'Tetraodon_nigroviridis', 'Oryza_sativa', 'Canis_familiaris', 'Latimeria_chalumnae', 'Geospiza_fortis', 'Musca_domestica']
    nameranking = np.zeros(len(protspecies))
    for s, spec in enumerate(species):
        if spec in speciesrank:
            nameranking[protspecies == s] = speciesrank.index(spec)
        else:
            nameranking[protspecies == s] = len(species)
    repclustnames = []
    bestclustmeasurement = []
    repdomains = []
    repdomainclass = []
    withinidentity = []
    for c, cl in enumerate(clusterassign):
        clmask = clusterid == c
        exnames = realprotnames[clmask][np.argsort(nameranking[clmask])]
        exdoms = np.copy(exnames)
        for e, en in enumerate(exnames):
            exnames[e], exdoms[e] = en.rsplit('__',1)
        exnames = exnames[:3]
        exdoms, edomnum = np.unique(exdoms,return_counts = True)
        exdoms = exdoms[np.argmax(edomnum)]
        exmeas, exmeasn = np.unique(closestmeas[clmask], return_counts = True)
        exmeas = exmeas[np.argmax(exmeasn)]
        repclustnames.append(','.join(exnames))
        bestclustmeasurement.append(exmeas)
        repdomains.append(exdoms)
        repdomainclass.append(exdoms.split('-')[0])
        if np.sum(clmask) > 1:
            withinidentity.append(np.mean(protsim[clmask][:,clmask][np.triu_indices(int(np.sum(clmask)))]))
        else:
            withinidentity.append(100)
    
    return repclustnames, repdomains, repdomainclass, withinidentity, bestclustmeasurement






















# Functions for visualization



def gainsandloss_at_branchpoints(unking, uniqueage, clusterage, cluster, evodistance, kings, shortspec, seqage, evolutiontail):
    
    xnames = []
    gainsplit = []
    evolution = []
    gains = []
    gainage = []
    gaininclude = []
    maintains = []
    losses = []
    lossspec = []
    colors = []
    kingcover = []
    speciscolors = -np.ones(len(clusterage))
    
    # characterize gains at losses for each time point
    gotking = np.zeros(len(unking))
    co = 0
    for ua, uniage in enumerate(uniqueage):
        if uniage > 1:
            ageclusts = np.where(clusterage == uniage)[0]
            cladespec = np.sum(cluster[ageclusts], axis = 0) > 0
            ## use evodistance to add other species that have less distance to cluster
            cladespec = np.sum(evodistance[cladespec] <= uniage, axis = 0)> 0
            unkings = np.unique(kings[cladespec])
            # only follow tree until unique kingdom is assigned
            if len(unkings) > 1 or gotking[unking==unkings[0]][0] == 0:
                # this if for visualization with imshow
                kcov = np.zeros(len(unking))
                kcov[np.isin(unking, unkings)] = 1
                kingcover.append(kcov)
                if co < 3:
                    rgbcolor = cm.tab20c_r(co+1)
                elif co < 23:
                    rgbcolor = cm.tab20b(co-3)
                else:
                    rgbcolor = cm.tab20(co-23)
                speciscolors[clusterage == uniage] = co
                colors.append(rgbcolor)
                co += 1
                xname = np.chararray(len(shortspec))
                xname[:] = '-'
                xname[cladespec] = shortspec[cladespec]
                xnames.append(str(uniage)+' MYa')
                gains.append(len(ageclusts))
                gainage.append(uniage)
                
                ages = np.array(gainage)
                numperage = np.zeros(len(ages))
                agess, agenum = np.unique(-seqage[ageclusts], return_counts= True)
                agess = -agess
                numperage[np.isin(ages,agess)] = agenum[np.isin(agess,ages)]
                if (agess < uniage).any():
                    numperage[-1] += np.sum(agenum[agess<uniage])
                gainsplit.append(numperage)
                
                evot, evotn = np.unique(evolutiontail[ageclusts], return_counts = True)
                evol = np.zeros(5)
                evol[evot] = evotn
                evolution.append(evol)
                
                gaininclude.append(cladespec)
                lossc = []
                loss = []
                maint = []
                # find last age that this cluster was part of
                compage = []
                for gi, ga in enumerate(gains):
                    if np.sum(gaininclude[gi]*cladespec) > 0:
                        compage.append(gi)
                for gi, ga in enumerate(gains):
                    if np.sum(gaininclude[gi]*cladespec) > 0:
                        gainclusts = cluster[clusterage == gainage[gi]] >=1
                        lostclusts = np.sum(np.sum(gainclusts * cladespec, axis = 1) == 0)
                        loss.append(lostclusts)
                        if len(compage)>1 and gi < len(gains)-1:
                            lossc.append(lostclusts - losses[compage[-2]][gi])
                        else:
                            lossc.append(0)
                        maint.append(ga-lostclusts)
                    else:
                        lossc.append(0)
                        loss.append(0)
                        maint.append(0)
                lossspec.append(lossc)
                losses.append(loss)
                maintains.append(maint)
                if len(unkings) == 1:
                    gotking[unking==unkings[0]] =1
    
    speciscolors[speciscolors == -1] = np.amax(speciscolors)+1
    # color for species specific
    colors.append(cm.BrBG(1.)) #(0.4,0.1,0.4,1.))
    return xnames, gainsplit, evolution, gains, gainage, gaininclude, maintains, losses, lossspec, colors, kingcover, speciscolors

def reduce_cluster_to_kingdom(cluster, speciscolors, unking, kings):
    kingcluster = np.zeros((len(cluster), len(unking)))
    for spco in np.unique(speciscolors)[:-1]:
        spcowhere = np.where(speciscolors == spco)[0]
        for c, csp in enumerate(cluster[spcowhere]):
            ckings = np.unique(kings[csp>=1])
            for ck in ckings:
                kingcount = csp[kings == ck]
                kingcluster[spcowhere[c], list(unking).index(ck)] = np.mean(kingcount[kingcount>0])
    return kingcluster
    
def plot_brachppoints(xnames, lossspec, gainsplit, evolution, colors, gains, kingcover, unking, outname, figsize = (4.5,2.5)):
    
    total = []
    for x, xname in enumerate(xnames):
        total.append(np.sum(lossspec[x]))
        total.append(np.sum(gainsplit[x]))
    
    for x, xname in enumerate(xnames):
        fig =plt.figure(figsize = figsize, dpi = 250)
        ax = fig.add_subplot(211)
        ax.set_position([0.35, 0.5, 0.3, 0.3])
        ax.set_title(xname)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_linewidth(1.4)
        ax.tick_params(labelright = False, labelleft = False, left = False, labelbottom = False, bottom = False, width = 1.4)
        ax.barh([2.5], np.cumsum(evolution[x])[::-1], color = cm.tab20([15, 8, 6, 10, 12]), lw = 0. , edgecolor = 'grey')
        ax.barh([1.5], np.cumsum(gainsplit[x])[::-1], color = colors[:x+1][::-1], lw = 0. , edgecolor = 'grey')
        lossx = np.cumsum(lossspec[x])[::-1]
        ax.barh([0], lossx, color = colors[:x+1][::-1])
        ax.text(np.amax(lossx)+1,0, str(-int(np.amax(lossx))), ha = 'left', va = 'center')
        ax.text(np.sum(gainsplit[x])+1,1.5, str(int(np.sum(gainsplit[x]))), ha = 'left', va = 'center')
        ax.set_ylim([-0.8, 3.3])
        ax.set_xlim([0, np.amax(total)])
        
        axp = fig.add_subplot(121)
        axp.set_position([0.2, 0.55, 0.12*2.5/4.5, 0.12])
        axp.set_axis_off()
        circle = plt.Circle((0.,0.), 1., facecolor = colors[x], ls = '-', lw = 0.75, edgecolor = 'k')
        axp.text(0.,0.,str(gains[x]), color = 'k', ha = 'center', va = 'center')
        axp.set_xlim([-1.1,1.1])
        axp.set_ylim([-1.1,1.1])
        axp.add_patch(circle)
        
        axc = fig.add_subplot(212)
        axc.set_position([0.3, 0.3, 0.5, 0.1])
        axc.tick_params(left = False, labelleft = False, bottom = False, which = 'both')
        axc.imshow([kingcover[x]], cmap = 'Greys', origin = 'lower', aspect = 'auto', vmin = 0, vmax = 1)
        axc.set_xticks(np.arange(len(unking)))
        axc.set_xticks(np.arange(len(unking))-0.5, minor = True)
        axc.grid(which = 'minor', color = 'white')
        axc.set_xticklabels(unking, rotation = 90, va = 'top', ha = 'center')
        
        fig.savefig(outname+xname.replace(' ', "_")+'.png', transparent =True, bbox_inches = 'tight', dpi = 200)
        print(outname+xname.replace(' ', "_")+'.png')
        plt.close()


if __name__ == '__main__':
    #python3 CRMG_clustering.py --analyze species_handcraft_measured4.nwk.txt -1 Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent.npz 0.2 ../JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz species_handcraft_measured4.nwkdistance_matrix.txt Orthogroups/species_handcraft_measured4_RRMKH_domains_concatenated_lnormminid.npz 70 30 --correct_domaincount --correct_length Orthogroups/species_handcraft_measured4_RRMKH_domains_concatenated.fasta 55 --allidlink
    #python3 CRMG_clustering.py --datasheet JPLEcluster_species_handcraft_measured4.nwk-1_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent_clusterscut0.2_idlink.txt species_handcraft_measured4.nwk.txt -1 species_handcraft_measured4.nwkdistance_matrix.txt --savefig --plot_branchpoints --evopercentage
    if '--analyze' in sys.argv:
        # list of species and their kingdoms
        species = np.genfromtxt(sys.argv[sys.argv.index('--analyze')+1], dtype = str)
        # column where to find the kinddom description
        kings = species[:, int(sys.argv[sys.argv.index('--analyze')+2])]
        species = species[:,0]
        # file names for latent representations
        latentrep = sys.argv[sys.argv.index('--analyze')+3]
        # cutoff value for clustering on latent distances
        confcut = float(sys.argv[sys.argv.index('--analyze')+4])
        
        # latent representation of measured RBPs
        latentmeas = np.load(sys.argv[sys.argv.index('--analyze')+5])
        measprots = latentmeas['names'].astype(str)
        measlatent = latentmeas['profiles']
        
        # distancce matrix between species
        evodistance = np.genfromtxt(sys.argv[sys.argv.index('--analyze')+6])
        evodistance = np.around(evodistance, 3)
        evospecies = open(sys.argv[sys.argv.index('--analyze')+6], 'r').readline().strip().replace("'",'').split()[1:]
        sort = []
        for s, spec in enumerate(evospecies):
            sort.append(list(species).index(spec))
        kings = kings[sort]
        species = species[sort]
        
        # sequence identity matrix for all RBPs used in the analysis
        protein_simfile = np.load(sys.argv[sys.argv.index('--analyze')+7])
        # cutoff identity that corrects for seperated RBPs (70%)
        
        psimcut = float(sys.argv[sys.argv.index('--analyze')+8])
        # cutoff distance below which RBP cannot be added to the cluster if not any other RBP has atleast this seqID to it
        plowcut = float(sys.argv[sys.argv.index('--analyze')+9])
        
        if '--outname' in sys.argv: 
            outname = sys.argv[sys.argv.index('--outname')+1]+os.path.splitext(os.path.split(sys.argv[sys.argv.index('--analyze')+1])[1])[0]+sys.argv[sys.argv.index('--analyze')+2]+'_'+os.path.splitext(sys.argv[sys.argv.index('--analyze')+3])[0]+'_clusterscut'+sys.argv[sys.argv.index('--analyze')+4]
        else:
            outname = 'JPLEcluster_'+os.path.splitext(os.path.split(sys.argv[sys.argv.index('--analyze')+1])[1])[0]+sys.argv[sys.argv.index('--analyze')+2]+'_'+os.path.splitext(sys.argv[sys.argv.index('--analyze')+3])[0]+'_clusterscut'+sys.argv[sys.argv.index('--analyze')+4]
        
        
        
        # Read in protein identities
        protsim = protein_simfile['identmat']
        protsimnames = protein_simfile['names'].astype(str)
        
        # generate arrays that determine which proteins can be compared to each other based on their domain composition
        # for example 1 rrm can be compared to 1 and 2 rrms but not 3
        maxdom = 5
        protsimnames, realprotnames, domaincounts, protdomains, protdomains2, remove = make_protdomain_arrays(protsimnames, maxdom=maxdom)
        # remove human pseudo proteins
        protsim = protsim[remove][:, remove]
        
        
        if '--correct_domaincount' in sys.argv:
            print('Correcting SeqID with domain count')
            # Generate matrix with nominator the minimum domain number of two RBPs and denominator the maximum of the two RBPs, set proteins that have x-1, x, x+1 domains to 1.
            domaincountmat = np.ones((len(domaincounts), len(domaincounts))) * domaincounts
            domaincountmat = np.amin(np.array([domaincountmat, domaincountmat.T]),axis = 0)/np.amax(np.array([domaincountmat, domaincountmat.T]),axis = 0)
            domaincountmat[np.dot(protdomains.T, protdomains2)>0] = 1.
            
            # correct sequence identity by minimum through maximum domain numbers
            protsim = protsim*domaincountmat

        
        # read in latent representations
        latent = []
        latentnames = []
        protspecies = []
        for s, spec in enumerate(species):
            slat = np.load(spec+'/'+latentrep)
            latent.append(slat['profiles'])
            latentnames.append(slat['names'].astype(str))
            protspecies.append(np.ones(len(latentnames[-1]), dtype = int)*s)
        
        
        # sort latent representation to protein similarity 
        latent = np.concatenate(latent, axis = 0)
        latentnames = np.concatenate(latentnames)
        protspecies = np.concatenate(protspecies)
        
        latsort = np.argsort(latentnames)[np.isin(np.sort(latentnames), protsimnames)]
        latentnames = latentnames[latsort]
        protspecies = protspecies[latsort]
        latent = latent[latsort]
        
        protsort = np.argsort(protsimnames)[np.isin(np.sort(protsimnames), latentnames)]
        protsimnames = protsimnames[protsort]
        realprotnames = realprotnames[protsort]
        protsim = protsim[protsort][:, protsort]
        protdomains = protdomains[:, protsort] >=1
        protdomains2 = protdomains2[:, protsort] >=1


        # remove proteins with short sequences, may be false detections from pfam scanning
        if '--correct_length' in sys.argv:
            lcorfile = sys.argv[sys.argv.index('--correct_length')+1]
            lcorlength = int(sys.argv[sys.argv.index('--correct_length')+2])
            lcorfile = open(lcorfile, 'r').readlines()
            lcrname = []
            lclen = []
            for l, line in enumerate(lcorfile):
                if line[0] == '>':
                    lcrname.append(line[1:].split('__')[0])
                else:
                    lclen.append(len(line.strip()))
            lclen = np.array(lclen)
            lcrname = np.array(lcrname)
            lcrname = lcrname[lclen >= lcorlength]
            sort = np.isin(latentnames, lcrname)
            
            print('Correct for minimum sequence length', lcorlength)
            print(len(latentnames))
            latentnames = latentnames[sort]
            protspecies = protspecies[sort]
            latent = latent[sort]
            protsimnames = protsimnames[sort]
            realprotnames = realprotnames[sort]
            protsim = protsim[sort][:, sort]
            protdomains = protdomains[:, sort]
            protdomains2 = protdomains2[:, sort]
            print(len(latentnames))

        # generate 3 matrices that mask only rbps that can be compared depending on their domain composition
        # has ones for two RBPs that have equal number of domains or one difference
        protdommat = np.dot(protdomains2.T, protdomains) > 0
        
        # arrays that have ones at position of the number of domains in the RBP
        protdomains = np.array([np.sum(protdomains[:maxdom], axis = 0) > 0, np.sum(protdomains[maxdom:], axis = 0) > 0])
        
        # is a matrix that has one if two rbps have exact same number of domains but doesn't matter which domiain type.
        protdommat2 = np.dot(protdomains.T, protdomains) > 0
        # loosest matrix allows for ones if dn1-1 = dn2+1
        protdommat3 = np.dot(protdomains2.T, protdomains2) > 0
        
        
        # compute latent distances to measured rbps
        latdist = cdist(measlatent, latent, 'cosine')
        closestdist = np.amin(latdist, axis = 0)
        closestmeas = measprots[np.argmin(latdist, axis = 0)]
        # compute latent distances between rbps
        flatdist = cdist(latent, latent, 'cosine')
        '''
        print(flatdist[latentnames == 'P103080_1.95d', latentnames == 'P086417_1.95d'])
        measdist = cdist(measlatent, measlatent, 'cosine')
        print(measdist[measprots == 'RNCMPT00127'][:, measprots == 'RNCMPT00419'])
        print(measdist[measprots == 'RNCMPT00127'][:, measprots == 'RNCMPT01680'])
        print(measdist[measprots == 'RNCMPT00058'][:, measprots == 'RNCMPT00419'])
        print(measdist[measprots == 'RNCMPT00058'][:, measprots == 'RNCMPT01680'])
        '''
        
        # compute homolog and ortholog matrix
        hmatrix, omatrix, tmatrix = homologmatrix(protsim, protspecies, testortholog = True, testtriangle= True, maskmat = protdommat)
        
        
        
        #perform agglomerative clustering with linkage
        # use protdommat which indicates if rbps with same domaintype can be assinged into same cluster
        clusterid = conserved_feature_group_clustering(hmatrix, flatdist,protspecies, confcut, maskmat = protdommat)
        clusterassign, unnum = np.unique(clusterid, return_counts = True)
        print('Unique cluster detected after evolution clustering', len(clusterassign), np.amax(clusterid))
        ccsizes, ccnsizes = np.unique(unnum, return_counts = True)
        print('Size distribution', '\n', ccsizes, '\n', ccnsizes)
                
        clusterid = combine_paralog_groups(clusterid, flatdist, protspecies, protdommat*(protsim >=30))
        clusterassign, unnum = np.unique(clusterid, return_counts = True)
        print('Unique cluster detected after paralog consideration', len(clusterassign), np.amax(clusterid))
        ccsizes, ccnsizes = np.unique(unnum, return_counts = True)
        print('Size distribution', '\n', ccsizes, '\n', ccnsizes)
        
        
        # make cluster matrix that contains the number of proteins per species
        # and assign cluster an age from the composition of species
        cluster, clusterage = assign_clusterage(clusterid, protspecies, evodistance)
        # need clusterage for fp correction
        
        # remove false positives from cluster if they harbour species outside the core of the cluster and have no support of being a one-to-one ortholog
        cluster, clusterage, clusterid = cleanfpcluster(tmatrix, evodistance, clusterage, cluster, clusterid, protspecies)        
        clusterassign, unnum = np.unique(clusterid, return_counts = True)
        print('Unique cluster after FP correction', len(clusterassign))
        ccsizes, ccnsizes = np.unique(unnum, return_counts = True)
        print('Size distribution', '\n', ccsizes, '\n', ccnsizes)
        
        # care for false negatives, which are all rbps that are 70% identical and also have same domain structure. 
        clusterid  = cleanfncluster(protsim, clusterid, psimcut, cluster)
        clusterassign, unnum = np.unique(clusterid, return_counts = True)
        print('Unique cluster after FN correction', len(clusterassign))
        ccsizes, ccnsizes = np.unique(unnum, return_counts = True)
        print('Size distribution', '\n', ccsizes, '\n', ccnsizes)
        
        # make cluster matrix that contains the number of proteins per species
        # and assign cluster an age from the composition of species
        cluster, clusterage = assign_clusterage(clusterid, protspecies, evodistance)
        
        
        # determine family relationship choosing from different methods
        if '--allidlink' in sys.argv:
            outname += '_idlink'
            family_relation, familyage, family, famevidence = determine_family(cluster, clusterage, evodistance, clusterid, protsim, trust= plowcut)
        elif '--homolink' in sys.argv:
            outname += '_homolink'
            family_relation, familyage, family, famevidence = determine_family(cluster, clusterage, evodistance, clusterid, (protsim+protsim*(tmatrix.astype(float)+omatrix.astype(float))/2.)/2., trust= None)
        else:
            family_relation, familyage, family, famevidence = determine_family(cluster, clusterage, evodistance, clusterid, protsim, trust= plowcut, maskmat = hmatrix)
        
        
        # Set all -1 clusters to their own number
        family[family == -1] = clusterassign[family == -1]
        
        age, agedist = np.unique(family_relation, return_counts = True)
        print('Final age distribution of RSSGs\n', age, '\n', agedist)
        
        # which clusters have a protein that is close to a measured RBP
        clustermeasbool = np.isin(clusterid, clusterid[closestdist < confcut])

        # generate per protein information matrices from cluster information
        clusterprot = np.zeros((len(clusterid), len(cluster[0])), dtype = np.int8)
        clusterages = np.zeros(len(clusterid))
        seqage = np.zeros(len(clusterid))
        evolutionorigin = np.zeros(len(clusterid), dtype = int)
        evolutionids = np.zeros(len(clusterid))
        evolutiontype = np.zeros(len(clusterid), dtype = np.int8)
        for c, cl in enumerate(clusterassign):
            clusterprot[clusterid == cl,:] = cluster[cl]
            clusterages[clusterid == cl] = clusterage[cl]
            seqage[clusterid == cl] = familyage[cl]
            evolutionorigin[clusterid == cl] = family[cl]
            evolutionids[clusterid == cl] = famevidence[cl]
            evolutiontype[clusterid == cl] = family_relation[cl]
        
        
        np.savetxt(outname+'.txt', np.append(np.array([latentnames, realprotnames, clusterid, clustermeasbool.astype(int), np.around(closestdist,3), closestmeas, np.sum(clusterprot>=1,axis = 1).astype(int), clusterages, seqage, evolutionorigin, evolutionids, evolutiontype, protspecies]).T, clusterprot, axis = 1), header = 'ProtId Protname Cluster Measured Closest_jple CloestRBR N_species Age_spec Age_seq Parent ParentSID Evolution SpeciesId; '+' '.join(species), fmt = '%s')
        print('SAVED protein specific information in', outname+'.txt')
        
        
        # Generate text file that is only for clusters and give one to three RBP names for each cluster: Names from certain species are preferred.
        
        # find cluster representatives
        # and best closest measurement for cluster
        repclustnames, repdomains, repdomainclass, withinidentity, bestclustmeasurement = generate_cluster_stats(clusterassign, clusterid, realprotnames, protspecies, species)

        # Generate a general name for the clade of clades
        assignedclade = assignclade(cluster)
        
        # Determine to how many groups this groups is a relative
        berelative, berelativenum = np.unique(family.astype(int), return_counts = True)
        isrelative = np.zeros(len(family), dtype = int)
        isrelative[berelative] = berelativenum
        
        # Determine if any RBP in the cluster is close to a measured RBP
        clustermeasbool = np.isin(np.arange(len(cluster)),np.unique(clusterid[clustermeasbool]))
        
        np.savetxt(outname+'_cluster.txt', np.append(np.array([clusterassign, clustermeasbool, np.array(repclustnames), np.array(repdomains), np.array(repdomainclass), np.array(bestclustmeasurement), np.sum(cluster,axis = 1).astype(int), np.around(withinidentity,0), clusterage, familyage, family.astype(int), np.around(famevidence,0), family_relation, np.sum(cluster>0,axis = 1), np.array(assignedclade), isrelative]).T, cluster, axis = 1), header = 'ClusterId;Measured;Proteinnames;Domainclass;Domaintype;RNAcmpt;Num_RBPs;IDwithin;Age_spec;Age_seq;RelativeRSSG;RelativeSID;Evolution;Num_species;Clade;Num_ancestor;'+';'.join(species), fmt = '%s', delimiter = ';')
        print('STATS saved as', outname+'_cluster.txt')
        
        sys.exit()
        
    elif '--datasheet' in sys.argv:
        # reads in protein specific datasheet to make figures
        
        #python3 CRMG_clustering.py --datasheet JPLEcluster_species_handcraft_measured4.nwk-1_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent_clusterscut0.2_idlink.txt species_handcraft_measured4.nwk.txt -1 species_handcraft_measured4.nwkdistance_matrix.txt --savefig --plot_branchpoints --evopercentage
        
        species = open(sys.argv[sys.argv.index('--datasheet')+1], 'r').readline().strip().split()[14:]
        kingfile = np.genfromtxt(sys.argv[sys.argv.index('--datasheet')+2], dtype = str)
        clade = int(sys.argv[sys.argv.index('--datasheet')+3])
        kings = []
        for s, spec in enumerate(species):
            kings.append(kingfile[list(kingfile[:,0]).index(spec),clade])
        kings = np.array(kings)
        
        evodistance = np.genfromtxt(sys.argv[sys.argv.index('--datasheet')+4])
        evodistance = np.around(evodistance, 3)
        evospecies = open(sys.argv[sys.argv.index('--datasheet')+4], 'r').readline().strip().replace("'",'').split()[1:]
        sort = []
        for s, spec in enumerate(species):
            sort.append(list(evospecies).index(spec))
        evodistance = evodistance[sort][:,sort]
        
        clusterprot = np.genfromtxt(sys.argv[sys.argv.index('--datasheet')+1], dtype = str)
        clusterprot = clusterprot[:,1:]
        clusterrealname = clusterprot[:,0]
        clusterid = clusterprot[:, 1].astype(int)
        clusterassign, clusterindex = np.unique(clusterid, return_index = True)
        clustermeasbool = clusterprot[clusterindex,2] == '1'
       
        if '--outname' in sys.argv:
            outname = sys.argv[sys.argv.index('--outname')+1]
        else:
            outname = os.path.splitext(sys.argv[sys.argv.index('--datasheet')+1])[0]
        cluster = clusterprot[:,12:].astype(int)[clusterindex]
        closestdist = clusterprot[:,3].astype(float)
        closestmeasured = clusterprot[:,4]
        clustermeasured = cluster[clustermeasbool]
        print('Measured', len(clustermeasured))
        clusterlen = clusterprot[:,5].astype(int)
        clusterage = clusterprot[:,6].astype(float)[clusterindex]
        seqage = clusterprot[:,7].astype(float)[clusterindex]
        protspecies = clusterprot[:,11].astype(int)
        evolutiontype = clusterprot[:,10].astype(int)
        evolutiontail = evolutiontype[clusterindex]
        
    else:
        print('Provide data sheet or analyze latent space')
        sys.exit()
    
    
    # Only generate plots for confident predictions
    if '--confident' in sys.argv:
        confident = float(sys.argv[sys.argv.index('--confident')+1])
        
        confidentclusters = np.unique(clusterid[closestdist <= confident])
        nonconfidentclusters = clusterassign[~np.isin(clusterassign, confidentclusters)]
        
        cid, cidd = np.unique(clusterid, return_index = True)
        clustlen = clusterlen[cidd]
        print('Confident Noncofident')
        print('Mean', np.mean(clustlen[confidentclusters]), np.mean(clustlen[nonconfidentclusters]))
        print('Number', len(confidentclusters), len(nonconfidentclusters))
        
        if '--nonconfident' in sys.argv:
            outname += '_nonconfidence'+str(confident)
            nonconfidentclusters, confidentclusters = confidentclusters, nonconfidentclusters
        else:
            outname += '_confidence'+str(confident)
        
    # Only generate plots for chosen evolutions
    if '--evolution' in sys.argv:
        chosenevos = np.array(sys.argv[sys.argv.index('--evolution')+1].split(','), dtype = int)
        
        confidentclusters = np.unique(clusterid[np.isin(evolutiontype, chosenevos)])
        nonconfidentclusters = clusterassign[~np.isin(clusterassign, confidentclusters)]
        
        cid, cidd = np.unique(clusterid, return_index = True)
        clustlen = clusterlen[cidd]
        print('Confident Noncofident')
        print('Mean', round(np.mean(clustlen[confidentclusters]),2), round(np.mean(clustlen[nonconfidentclusters]),2))
        print('Number', len(confidentclusters), len(nonconfidentclusters))
        
        if '--nonconfident' in sys.argv:
            outname += '_nonconfidence'+str(confident)
            nonconfidentclusters, confidentclusters = confidentclusters, nonconfidentclusters
        else:
            outname += '_evolution'+'-'.join(chosenevos.astype(str))

    
    if '--domainspecific' in sys.argv:
        domainspec = sys.argv[sys.argv.index('--domainspecific')+1]
        realdomain = []
        for c, clrna in enumerate(clusterrealname):
            realdomain.append(clrna.split('__')[-1])
        realdomain = np.array(realdomain)
        confidentclusters = np.unique(clusterid[realdomain == domainspec])
        nonconfidentclusters = clusterassign[~np.isin(clusterassign, confidentclusters)]
        cid, cidd = np.unique(clusterid, return_index = True)
        clustlen = clusterlen[cidd]
        print(domainspec)
        print('Mean', np.mean(clustlen[confidentclusters]), np.mean(clustlen[nonconfidentclusters]))
        print('Number', len(confidentclusters), len(nonconfidentclusters))
        outname += '_'+domainspec
    
    
    if '--domainspecific' in sys.argv or '--confident' in sys.argv or '--evolution' in sys.argv:
        cluster = cluster[np.isin(clusterassign,confidentclusters)]
        clustermeasbool = clustermeasbool[np.isin(clusterassign,confidentclusters)]
        #clustercenter = clustercenter[np.isin(clusterid,confidentclusters)]
        closestdist = closestdist[np.isin(clusterid,confidentclusters)]
        clusterlen = clusterlen[np.isin(clusterid,confidentclusters)]
        clusterage = clusterage[np.isin(clusterassign,confidentclusters)]
        
        clusterrealname = clusterrealname[np.isin(clusterid,confidentclusters)]
        evolutiontype = evolutiontype[np.isin(clusterid,confidentclusters)]
        evolutiontail = evolutiontail[np.isin(clusterassign,confidentclusters)]
        seqage = seqage[np.isin(clusterassign,confidentclusters)]
        
        clusterassign = clusterassign[np.isin(clusterassign,confidentclusters)]
        clusterid = clusterid[np.isin(clusterid,confidentclusters)]
        for c, clusa in enumerate(clusterassign):
            clusterid[clusterid == clusa] = c
        clusterassign = np.unique(clusterid)
        
    
    speciclusters = np.sum(cluster>0, axis = 0)
    speciclustersmeasured = np.sum(clustermeasured, axis = 0)
    
    # generate unique set of kingdoms
    unking = []
    for kun in kings:
        if kun not in unking:
            unking.append(kun)
    unking = np.array(unking)
    
    # unique ages of clusters
    uniqueage, ageid = np.unique(clusterage, return_index = True)
    uniqueage, ageid = uniqueage[::-1], ageid[::-1]
    
    # short name for species
    shortspec = []
    for s, spec in enumerate(species):
        shortspec.append(spec[0])
    shortspec = np.array(shortspec)
    
    # Compute gains and losses and their respected colours in the plots
    xnames, gainsplit, evolution, gains, gainage, gaininclude, maintains, losses, lossspec, colors, kingcover, speciscolors = gainsandloss_at_branchpoints(unking, uniqueage, clusterage, cluster, evodistance, kings, shortspec, seqage, evolutiontail)
    
    # reduce cluster assignment to kingdoms
    kingcluster = reduce_cluster_to_kingdom(cluster, speciscolors, unking, kings)
    kingclustermeasured = kingcluster[clustermeasbool]


    # plot statistics for brachpoints
    if '--plot_branchpoints' in sys.argv:
        plot_brachppoints(xnames, lossspec, gainsplit, evolution, colors, gains, kingcover, unking, outname, figsize = (4.5,2.5))
        
        
    #### MADE MORE readable until here! Need to continue
    
    
    speciesorigin = []
    sequence_age = []
    seqagecolors = []
    
    speciesevolution = []
    speciesbars = []
    speciesbarslost = []
    species_colors = []
    species_colorslost = []
    species_age = []
    
    # compute all different color assignments and how often
    totalcluster, totalind, totalnum = np.unique(speciscolors, return_counts = True, return_index = True)
    totalages = clusterage[totalind.astype(int)]
    totalages[-1] = 0.
    
    
    for s, spec in enumerate(species):
        # find all different types of evolution colors
        speiclss = speciscolors[cluster[:,s] >= 1]
        # determine unique evolution colors in species set
        uncolor, unind, unnum = np.unique(speiclss, return_index = True, return_counts = True)
        # compuate all clusterages that are in the youngest color categorie
        addcolor, addcolornum = np.unique(-clusterage[(cluster[:,s] >=1)*(speciscolors == uncolor[-1])], return_counts = True)
        # collect colors for lost specificities, don't have the color of single ones
        species_colorslost.append(np.array(colors)[uncolor[:-1].astype(int)][::-1])
        # change the darkness of color for youngest class, to distinguish between different subgroups that come after the clade name
        # split last group up into different categories
        morenum = np.append(unnum[:-1], addcolornum)
        #choose last color
        morecolor = np.array(colors)[uncolor.astype(int)][:-1]
        # for all ages in youngest class do
        for a, adc in enumerate(addcolor):
            changecolor = np.array(colors[-1])
            changecolor[-1] = 0.7 + 0.3*-adc/max(-addcolor[0],1.) #(len(addcolor)-1-a)/max(1,len(addcolor)-1)
            morecolor = np.append(morecolor, [changecolor], axis = 0)
        
        speciesbars.append(np.cumsum(morenum)[::-1])
        species_colors.append(morecolor[::-1])
        species_age.append(clusterage[cluster[:,s] >= 1][unind])
        
        #print(speciesbars[-1])
        
        seqagenum = np.zeros(len(totalcluster))
        useqage, useqagenum = np.unique(-seqage[cluster[:,s]>=1], return_counts = True)
        useqage = -useqage
        if np.sum(np.isin(useqage, totalages)) != len(useqage):
            u0 = np.sum(useqagenum[~np.isin(useqage, totalages)])
            useqagenum = useqagenum[np.isin(useqage, totalages)]
            useqage = useqage[np.isin(useqage, totalages)]
            if 0 in useqage:
                useqagenum[useqage == 0] += u0
            else:
                useqage = np.append( useqage, [0])
                useqagenum = np.append(useqagenum, [u0])
            
        seqagenum[np.isin(totalages, useqage)] = useqagenum
        sequence_age.append(np.cumsum(seqagenum[seqagenum>0])[::-1])
        seqagecolors.append(np.array(colors)[seqagenum>0][::-1])
        
        potages = list(np.unique(evodistance[s]))
        #print(spec, potages)
        specorn = np.zeros((len(potages), len(potages)))
        seqtocluster, staing, seqtoclustnum = np.unique(np.column_stack((seqage[cluster[:,s]>=1], clusterage[cluster[:,s]>=1])), return_counts = True, return_index= True, axis = 0)
        for tc, stc in enumerate(seqtocluster): 
            #print(stc, seqtoclustnum[tc])
            specorn[potages.index(stc[0]), potages.index(stc[1])] = seqtoclustnum[tc]
        speciesorigin.append(specorn)
        
        unnumspecies = totalnum[np.isin(totalcluster, uncolor)]-unnum 
        speciesbarslost.append(np.cumsum(unnumspecies[:-1])[::-1])
        sevot, sevotn = np.unique(evolutiontail[cluster[:, s]>=1], return_counts = True)
        evoln = np.zeros(5)
        evoln[sevot] = sevotn
        evoln[1], evoln[2] = evoln[2], evoln[1]
        speciesevolution.append(np.cumsum(evoln)[::-1])
        #print(speciesevolution[-1])
    
    

    # plot of retained number of specificities per species 
    if '--horizontalbar' in sys.argv:
        figs = plt.figure(figsize = (len(species)*0.28, 4), dpi = 100)
    else:
        figs = plt.figure(figsize = (4.5, len(species)*0.28), dpi = 100)
    axs = figs.add_subplot(111)
    for s, spec in enumerate(species):
        if '--horizontalbar' in sys.argv:
            axs.bar(np.ones(len(speciesbars[s]))*s, speciesbars[s], color = species_colors[s], linewidth = 0., edgecolor = 'k')
        else:
            axs.barh(np.ones(len(speciesbars[s]))*-s, speciesbars[s], color = species_colors[s], linewidth = 0., edgecolor = 'k')
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    if '--horizontalbar' in sys.argv:
        axs.set_xticks(np.arange(len(species)))
        axs.set_xticklabels(species, rotation = 90)
    else:
        axs.set_yticks(-np.arange(len(species)))
        axs.set_yticklabels(species, rotation = 0)
        axs.set_ylim([-len(species), 0.5])
        axs.set_xticks(np.arange(0, np.amax(np.concatenate(speciesbars)), 20, dtype = int))
        axs.grid(axis = 'x')
    axs.set_xlabel('Number specificity clusters')
    
    axslim = axs.get_xlim()
    
    # plot of lost number of specificities
    if '--horizontalbar' in sys.argv:
        figs2 = plt.figure(figsize = (len(species)*0.28, 4), dpi = 100)
    else:
        figs2 = plt.figure(figsize = (4.5, len(species)*0.28), dpi = 100)
    axs2 = figs2.add_subplot(111)
    for s, spec in enumerate(species):
        if '--horizontalbar' in sys.argv:
            axs2.bar(np.ones(len(speciesbarslost[s]))*s, speciesbarslost[s], color = species_colorslost[s], linewidth = 0., edgecolor = 'k')
        else:
            #print(speciesbarslost[s])
            axs2.barh(np.ones(len(speciesbarslost[s]))*-s, speciesbarslost[s], color = species_colorslost[s], linewidth = 0., edgecolor = 'k')
    axs2.spines['top'].set_visible(False)
    axs2.spines['right'].set_visible(False)
    if '--horizontalbar' in sys.argv:
        axs2.set_xticks(np.arange(len(species)))
        axs2.set_xticklabels(species, rotation = 90)
    else:
        axs2.set_yticks(-np.arange(len(species)))
        axs2.set_yticklabels(species, rotation = 0)
        axs2.set_ylim([-len(species), .5])
        axs2.set_xticks(np.arange(0, np.amax(np.concatenate(speciesbars)), 20, dtype = int))
        axs2.grid(axis = 'x')
    axs2.set_xlim(axslim)
    axs2.set_xlabel('Changed specificity clusters')
    
    # plot of evolutionary relationships of retained specificities
    if '--horizontalbar' in sys.argv:
        figs3 = plt.figure(figsize = (len(species)*0.28, 4), dpi = 100)
    else:
        if '--evopercentage' in sys.argv or '--evobinomial' in sys.argv:
            figs3 = plt.figure(figsize = (3., len(species)*0.28), dpi = 100)
        else:
            figs3 = plt.figure(figsize = (4.5, len(species)*0.28), dpi = 100)
            
    axs3 = figs3.add_subplot(111)
    for s, spec in enumerate(species):
        if '--evopercentage' in sys.argv:
            #print(speciesevolution[s])
            speciesevolution[s] = speciesevolution[s] - speciesevolution[s][-1]
            speciesevolution[s][speciesevolution[s]<0] = 0
            speciesevolution[s] = speciesevolution[s]*100./max(1,np.amax(speciesevolution[s]))
        elif '--evobinomial' in sys.argv:
            speciesevolution[s] = speciesevolution[s] - speciesevolution[s][-1]
            nsamples = speciesevolution[s][1]
            speciesevolution[s][speciesevolution[s]<0] = 0
            speciesevolution[s] = speciesevolution[s]*100./max(1,speciesevolution[s][1])
            speciesevolution[s][-3] = np.sqrt(((speciesevolution[s][-2])/100.*(100.-speciesevolution[s][-2])/100.)/nsamples)*100.
            
        if '--horizontalbar' in sys.argv:
            axs3.bar(np.ones(len(speciesevolution[s]))*s, speciesevolution[s], color = cm.tab20([15, 9,10,6,12]), linewidth = 0., edgecolor = 'k')
        else:
            if '--evobinomial' in sys.argv:
                axs3.barh(np.ones(len(speciesevolution[s]))*-s, speciesevolution[s][-2], xerr = speciesevolution[s][-3], color =cm.tab20(10), linewidth = 0., edgecolor = 'k', )
            else:
                axs3.barh(np.ones(len(speciesevolution[s]))*-s, speciesevolution[s], color =cm.tab20([15, 9,10,6,12]), linewidth = 0., edgecolor = 'k')
    axs3.spines['top'].set_visible(False)
    axs3.spines['right'].set_visible(False)
    if '--horizontalbar' in sys.argv:
        axs3.set_xticks(np.arange(len(species)))
        axs3.set_xticklabels(species, rotation = 90)
    else:
        axs3.set_yticks(-np.arange(len(species)))
        axs3.set_yticklabels(species, rotation = 0)
        axs3.set_ylim([-len(species), .5])
        if '--evopercentage' in sys.argv or '--evobinomial' in sys.argv:
            axs3.set_xticks(np.arange(0, 120, 20, dtype = int))
        else:
            axs3.set_xticks(np.arange(0, np.amax(np.concatenate(speciesevolution)), 20, dtype = int))
        axs3.grid(axis = 'x')
    
    if '--evopercentage' in sys.argv:
        axs3.set_xlabel('Percentage evolutionary relationship')
    elif '--evobinomial' in sys.argv:
        axs3.set_xlabel('Percentage duplication')
    else:
        axs3.set_xlabel('Number specificity clusters')
    
    
    
    # sequence ages
    if '--horizontalbar' in sys.argv:
        figs4 = plt.figure(figsize = (len(species)*0.28, 4), dpi = 100)
    else:
        figs4 = plt.figure(figsize = (4.5, len(species)*0.28), dpi = 100)
            
    axs4 = figs4.add_subplot(111)
    for s, spec in enumerate(species):
        if '--horizontalbar' in sys.argv:
            axs4.bar(np.ones(len(sequence_age[s]))*s, sequence_age[s], color = seqagecolors[s], linewidth = 0., edgecolor = 'k')
        else:
            axs4.barh(np.ones(len(sequence_age[s]))*-s, sequence_age[s], color =seqagecolors[s], linewidth = 0., edgecolor = 'k')
    axs4.spines['top'].set_visible(False)
    axs4.spines['right'].set_visible(False)
    if '--horizontalbar' in sys.argv:
        axs4.set_xticks(np.arange(len(species)))
        axs4.set_xticklabels(species, rotation = 90)
    else:
        axs4.set_yticks(-np.arange(len(species)))
        axs4.set_yticklabels(species, rotation = 0)
        axs4.set_ylim([-len(species), .5])
        axs4.set_xticks(np.arange(0, np.amax(np.concatenate(sequence_age)), 20, dtype = int))
        axs4.grid(axis = 'x')
    axs4.set_xlabel('Number specificity clusters')
    
    
    
    ########################################################
    #### NEW ########
    # generate the matrix
    gotking = np.zeros(len(unking))
    keep = []
    for u, uniage in enumerate(uniqueage):
        ageclusts = np.where(clusterage == uniage)[0]
        cladespec = np.sum(cluster[ageclusts], axis = 0) > 0
        ## use evodistance to add other species that have less distance to cluster
        cladespec = np.sum(evodistance[cladespec] <= uniage, axis = 0)> 0
        unkings = np.unique(kings[cladespec])
        # only follow tree until unique kingdom is assigned
        if len(unkings) > 1 or gotking[unking==unkings[0]][0] == 0:
            keep.append(True)
            if len(unkings) == 1:
                gotking[unking==unkings[0]] =1
        else:
            keep.append(False)
    
    uniage = uniqueage[keep][:-1]
    
    ageevmatrix = np.zeros((len(uniage), len(uniage)), dtype = int)
    for u, wniage in enumerate(uniage):
        for v, vniage in enumerate(uniage):
            ageevmatrix[v,u] = int(np.sum((seqage == wniage) * (speciscolors == v)))    
    
    for u in range(len(uniage)):
        if ageevmatrix[u,u] == 0:
            ageevmatrix[u,u] = -np.amax(ageevmatrix)
    
    
    cladenames = [bestmatchedcladename(unking[kingcover[c] > 0]) for c in range(len(kingcover))] 
    clna, clnn = np.unique(cladenames, return_counts = True)
    clad = np.ones(len(clna))
    clad[clnn > 1] = 0
    for c, cln in enumerate(cladenames):
        if clad[clna == cln] < clnn[clna == cln]:
            cladenames[c] += str(int(clad[clna == cln]))
            clad[clna == cln]+= 1
    
    soort = [cladenames.index(i) for i in ['LECA', 'LECA-some', 'Fungi+Metazoan+Plant', 'Green plants', 'Land plants0', 'Land plants1', 'Flowering_plant', 'Fungi+Metazoan', 'Metazoan', 'Vertebrates', 'Land-vertebrates0',  'Land-vertebrates1',  'Land-vertebrates2',  'Mammals0', 'Mammals1', 'Mammals2', 'Bird',  'Fish', 'Protostomia', 'Worm', 'Insects0','Insects1', 'Insect', 'Fly',  'Protists', 'Alveolata', 'Stramenopiles', 'Heterokonts', 'Plant_pathogen', 'Plasmodium']]
    ageevmatrix = ageevmatrix[soort][:, soort]
    cladenames = np.array(cladenames)[soort]
    uniage = uniage[soort]
    print(kingcover)
    kingcover = np.array(kingcover)[soort]
    
    # sequence ages versus specificity ages
    asi, asigap = 4, 1.5
    xunits, yunits = len(uniage)+len(kingcover[0])+asi, len(uniage)+asi
    oneunitx, oneunity = 0.8/xunits, 0.8/yunits
    figmatev = plt.figure(figsize = (xunits*0.3, yunits*0.3), dpi = 50)
    
    axmt = figmatev.add_subplot(131)
    axmt.set_position([0.1+(len(kingcover[0])+asi)*oneunitx, 0.1+asi*oneunity, len(uniage)*oneunitx ,0.8-asi*oneunity])
    axmt.spines['top'].set_visible(False)
    axmt.spines['left'].set_visible(False)
    axmt.spines['right'].set_visible(False)
    axmt.spines['bottom'].set_visible(False)
    axmt.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    
    axmt.imshow(np.sign(ageevmatrix) * np.log(1+np.absolute(ageevmatrix)), origin = 'upper', aspect = 'auto', cmap = 'BrBG')
    
    for i in range(len(ageevmatrix)):
        for j in range(len(ageevmatrix[0])):
            if ageevmatrix[i,j] > 0:
                axmt.text(j,i,str(ageevmatrix[i,j]), va = 'center', ha = 'center')
    
            
    
    
    axmt.set_xticks(np.arange(len(uniage)-1)+0.5, minor =True)
    axmt.set_yticks(np.arange(len(uniage))+0.5, minor =True)
    axmt.grid(which = 'minor')
    
    
    
    #cladenames = [bestmatchedcladename(unking[kingcover[c] > 0]) for c in range(len(kingcover))]
    axc = figmatev.add_subplot(132)
    axc.set_position([0.1, 0.1+asi*oneunity, len(kingcover[0])*oneunitx ,0.8-asi*oneunity])
    axc.tick_params(left = False, labelleft = True, bottom = False, which = 'both')
    axc.imshow(kingcover, cmap = 'Greys', origin = 'upper', aspect = 'auto', vmin = 0, vmax = 1)
    axc.set_xticks(np.arange(len(unking)))
    axc.set_xticks(np.arange(len(unking))-0.5, minor = True)
    axc.set_yticks(np.arange(len(kingcover))-0.5, minor = True)
    axc.set_yticks(np.arange(len(kingcover)))
    axc.set_yticklabels(cladenames)
    axc.grid(which = 'minor', color = 'grey', axis = 'x')
    axc.grid(which = 'minor', color = 'black', axis = 'y')
    axc.set_xticklabels(unking, rotation = 90, va = 'top', ha = 'center')
    
    # generate own colormap 
    lcolors = ListedColormap(colors)
     
    axd = figmatev.add_subplot(133)
    axd.set_position([0.1+(len(kingcover[0]) +asigap)*oneunitx, 0.1+asi*oneunity, (asi-asigap)*oneunitx ,0.8-asi*oneunity])
    axd.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False, which = 'both')
    axd.imshow(np.arange(len(uniage)).reshape(-1,1), cmap = lcolors, origin = 'upper', aspect = 'auto', vmin = 0, vmax = len(uniage))
    for i in range(len(uniage)):
        axd.text(0,i,str(round(uniage[i],1)), va = 'center', ha = 'center')
    axd.set_xlabel('MYA')
    axd.set_ylabel('RNA specificity age')
    
    axe = figmatev.add_subplot(144)
    axe.set_position([0.1+(len(kingcover[0])+asi)*oneunitx, 0.1+asigap*oneunity, len(uniage)*oneunitx ,(asi-asigap)*oneunity])
    axe.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False, which = 'both')
    axe.imshow(np.arange(len(uniage)).reshape(1,-1), cmap = lcolors, origin = 'lower', aspect = 'auto', vmin = 0, vmax = len(uniage))
    for i in range(len(uniage)):
        axe.text(i,0, str(round(uniage[i],1)), va = 'center', ha = 'center', rotation = 90)
    axe.set_xlabel('Sequence age')
    '''
    axf = figmatev.add_subplot(144)
    axf.set_position([0.9 + 0.5*oneunitx, 0.1+asi*oneunity,0.1,0.8-asi*oneunity])
    axf.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False, which = 'both')
    axf.spines['top'].set_visible(False)
    axf.spines['right'].set_visible(False)
    
    ageevmatrix[ageevmatrix <0] = 0
    cumulated = np.cumsum(ageevmatrix, axis = 1)/np.sum(ageevmatrix, axis = 1)[:,None]
    print(cumulated)
    cumulated = cumulated[:,::-1] * (ageevmatrix[:,::-1] > 0).astype(int)
    print(cumulated)
    
    for c in range(len(cumulated)):
        cumask = cumulated[c] > 0
        print(np.where(cumask)[0], len(colors))
        print(cumulated[c][cumask])
        axf.barh(-np.ones(int(np.sum(cumask)))*c, cumulated[c][cumask], color = np.array(colors)[::-1][1:][cumask], linewidth = 0., edgecolor = 'k')
    axf.set_ylim([-len(cumulated)+0.5, 0.5])
    axf.set_xlabel('Percentage')
    '''
    
    ##########################################################
    
    
    
    # need to work with maintains, gains, loss, gainage, gaininclude
    speciesbarslost = []
    speciesbars = []
    species_age = []
    totalcluster, totalnum = np.unique(-clusterage, return_counts = True)
    for s, spec in enumerate(species):
        unage = np.unique(-evodistance[s])
        specbar, speclos, specage = [], [], []
        for ua, uniage in enumerate(unage):
            cladespec = evodistance[s] <= -uniage
            gainclusts = cluster[np.isin(clusterage, -unage[:ua+1])] >=1
            clusts = np.sum(np.sum(gainclusts * cladespec, axis = 1) >= 1)
            specbar.append(clusts)
            specage.append(uniage)
            speclos.append(np.sum(totalnum[np.isin(totalcluster,unage[:ua+1])])-clusts)
            
        speciesbars.append(np.array(specbar))
        speciesbarslost.append(np.array(speclos))
        species_age.append(-np.array(specage))
        

    
    cladecolor = ['skyblue', 'blue', 'silver', 'indianred', 'sienna', 'green', 'limegreen', 'orange', 'grey', 'seagreen', 'purple']
    chosenspecies = ['Danio_rerio', 'Gallus_gallus', 'Homo_sapiens', 'Drosophila_melanogaster', 'Caenorhabditis_elegans', 'Arabidopsis_thaliana','Selaginella_moellendorffii']
    if '--separate_timeplot' not in sys.argv:
        fige = plt.figure(figsize = (6, 5.), dpi = 200)
        axe = fige.add_subplot(111)
        for q, spec in enumerate(chosenspecies):
            s = list(species).index(spec)
            axe.plot(-species_age[s], speciesbars[s], color = cladecolor[q], marker = 'o', markersize = 3., label = spec)
        axe.legend(prop={'size':5.5})
        axe.spines['top'].set_visible(False)
        axe.spines['right'].set_visible(False)
        axe.set_ylabel('Number specificities')
        #axe.set_xscale('symlog')
        axe.set_xticks([-1768, -1500, -1200, -900, -600, -300, 0])
        axe.set_xticklabels(['LECA', '1500 MYa', '1200 MYa', '900 MYa', '600 MYa', '300 MYa', '0 MYa'], rotation = 60)
        
        axe.plot([-np.amax(np.concatenate(species_age)), 0], [0,0], 'k', ls = '-')
        for q, spec in enumerate(chosenspecies):
            s = list(species).index(spec)
            print(spec,-speciesbarslost[s][:-1] )
            axe.plot(-species_age[s][:-1], -speciesbarslost[s][:-1], color = cladecolor[q], marker = 'v', markersize = 3., ls= '--', label = spec)
        print('Saved as', outname+'_species_clusterevolution.svg')
        fige.savefig(outname+'_species_clusterevolution.svg', bbox_inches = 'tight', dpi = 300)
        

    if '--separate_timeplot' not in sys.argv:
        speciesbarslost = []
        speciesbars = []
        species_age = []
        totalcluster, totalnum = np.unique(-clusterage, return_counts = True)
        for s, spec in enumerate(species):
            unage, unnum = np.unique(-clusterage[cluster[:,s]>=1], return_counts = True)
            speciesbars.append(np.cumsum(unnum))
            species_age.append(unage)
            unnumlost = totalnum[np.isin(totalcluster, unage)] - unnum
            speciesbarslost.append(np.cumsum(unnumlost[:-1]))
        
        fige = plt.figure(figsize = (5, 5), dpi = 200)
        axe = fige.add_subplot(111)
        for q, spec in enumerate(chosenspecies):
            s = list(species).index(spec)
            axe.plot(species_age[s], speciesbars[s], color = cladecolor[q], marker = 'o', markersize = 3., label = spec)
        axe.legend(prop={'size':5.5})
        axe.spines['top'].set_visible(False)
        axe.spines['right'].set_visible(False)
        axe.set_ylabel('Number specificities')
        #axe.set_xscale('symlog')
        axe.set_xticks([-1768, -1500, -1200, -900, -600, -300, 0])
        axe.set_xticklabels(['LECA', '1500 MYa', '1200 MYa', '900 MYa', '600 MYa', '300 MYa', '0 MYa'], rotation = 60)
        
        fige.savefig(outname+'_specificitytiming.svg', bbox_inches = 'tight', dpi = 300)
        
        
    else:
        for q, spec in enumerate(chosenspecies):
            fige = plt.figure(figsize = (4, 5.5), dpi = 200)
            axe = fige.add_subplot(111)    
            s = list(species).index(spec)
            axe.plot(-species_age[s], speciesbars[s], color = cladecolor[q], marker = 'o', label = spec)
            axe.legend(prop={'size':5.5})
            axe.spines['top'].set_visible(False)
            axe.spines['right'].set_visible(False)
            axe.set_xlabel('Mya')
            axe.set_ylabel('Number specificities')
            axe.set_ylim([0, np.amax(np.concatenate(speciesbars))+4])
            axe.set_xticks([-1768, -1500, -1200, -900, -600, -300, 0])
            axe.set_xticklabels(['LECA', '1500 MYa', '1200 MYa', '900 MYa', '600 MYa', '300 MYa', '0 MYa'], rotation = 60)
            fige.savefig(outname+'_species_clusterevolution'+spec+'.svg', bbox_inches = 'tight', dpi = 200)
            
   
    
    clusteroverlaps, clustid, barplot = np.unique(kingcluster>0, axis = 0, return_index = True, return_counts = True)
    clustage = clusterage[clustid]
    clusterparalogs = np.zeros(np.shape(clusteroverlaps))
    for co, clof in enumerate(clusteroverlaps):
        coset = []
        for k, kcu in enumerate(kingcluster >0):
            if np.array_equal(kcu, clof):
                coset.append(k)
        clusterparalogs[co] = np.mean(kingcluster[coset], axis = 0)
    
    
    barmeas = np.zeros(len(barplot))
    for c, clus in enumerate(kingclustermeasured>0):
        clus = clus >=1
        found = False
        for cv, clvo in enumerate(clusteroverlaps):
            if np.array_equal(clvo, clus):
                found = True
                break
        if found:
            barmeas[cv] += 1
        else:
            print("Not found")
            sys.exit()
    
    
    
    
    barmask = np.sum(clusteroverlaps, axis = 1) > 0
    clusterparalogs, clusteroverlaps, clustage, barmeas, barplot = clusterparalogs[barmask], clusteroverlaps[barmask], clustage[barmask], barmeas[barmask], barplot[barmask]
    
    barmask = np.argsort(-clustage+np.sum(clusteroverlaps,axis = 1))
    clusterparalogs, clusteroverlaps, clustage, barmeas, barplot = clusterparalogs[barmask], clusteroverlaps[barmask], clustage[barmask], barmeas[barmask], barplot[barmask]
    
    if '--filterrare' in sys.argv:
        barmask = (barplot>2) | ((barmeas>=1)&(barplot == 2))
        clusterparalogs, clusteroverlaps, clustage, barmeas, barplot = clusterparalogs[barmask], clusteroverlaps[barmask], clustage[barmask], barmeas[barmask], barplot[barmask]
    
    if '--filter3' in sys.argv:
        barmask = (barplot>3) | ((barmeas>=1)&(barplot == 3))
        clusterparalogs, clusteroverlaps, clustage, barmeas, barplot = clusterparalogs[barmask], clusteroverlaps[barmask], clustage[barmask], barmeas[barmask], barplot[barmask]
    
    if '--filterage' in sys.argv:
        barmask = clustage > 60
        clusterparalogs, clusteroverlaps, clustage, barmeas, barplot = clusterparalogs[barmask], clusteroverlaps[barmask], clustage[barmask], barmeas[barmask], barplot[barmask]
        
    
    if '--motif_evolution' in sys.argv:
        # get measured names that changed specificity to generate
        # motif matrix with confident motifs
        # seqid matrix
            # only for the ones that have measured ones
        # jple e-dist matrix
        for c, cla in enumerate(np.unique(clustage)):
            cltype = np.where(clustage == cla)[0]
            print(cltype)
            if len(cltype) > 1:
                cltype = cltype[:-1]
                for clt in cltype:
                    for k, kingc in enumerate(kingcluster>0):
                        if np.array_equal(kingc, clusteroverlaps[clt]):
                            print(k, clustercenter[clusterid == k][0], cla, ''.join(kingc.astype(int).astype(str)))
    if '--print_clusters' in sys.argv:
        speciesrank = ['Homo_sapiens', 'Mus_musculus', 'Macaca_mulatta', 'Gallus_gallus', 'Danio_rerio', 'Xenopus_tropicalis', 'Saccharomyces_cerevisiae', 'Schizosaccharomyces_pombe', 'Drosophila_melanogaster', 'Drosophila_ananassae', 'Caenorhabditis_elegans', 'Caenorhabditis_briggsae', 'Arabidopsis_thaliana', 'Arabidopsis_lyrata', 'Zea_mays', 'Tetraodon_nigroviridis', 'Oryza_sativa', 'Canis_familiaris', 'Latimeria_chalumnae', 'Geospiza_fortis', 'Musca_domestica']
        speciesdic = {}
        for s, spec in enumerate(species):
            if spec in speciesrank:
                speciesdic[s] = speciesrank.index(spec)
            else:
                speciesdic[s] = len(species)
        clustername = []
        for c, cla in enumerate(clusterassign):
            clnames = clusterrealname[clusterid == cla]
            ranking = [speciesdic[ps] for ps in protspecies[clusterid == cla]]
            clustername.append(clnames[np.argsort(ranking)][:3])
        
        
        print('clusterID', '#RBR', 'REP_RNCMPT', 'Ancestor', ','.join( unking), 'top10Names')
        for c, cla in enumerate(np.unique(clustage)):
            cltype = np.where(clustage == cla)[0]
            for clt in cltype:
                for k, kingc in enumerate(kingcluster):
                    if np.array_equal(kingc>0, clusteroverlaps[clt]):
                        bestrep = np.argmin(closestdist[clusterid == k])
                        #closerep = clusterrealname[np.where(clusterid == k)[0][np.argsort(closestdist[clusterid == k])]][:10]
                        print(k, int(np.sum(clusterid == k)), closestmeasured[clusterid == k][bestrep], cla, ','.join(np.around(kingc).astype(int).astype(str)), ','.join(clustername[k]))
        
       
    fig = plt.figure(figsize = (len(clusteroverlaps)*0.25, 6), dpi = 200)
    ax = fig.add_subplot(121)
    ax.set_position([0.1,0.5, 0.8,0.4])
    ax.bar(np.arange(len(barplot)), barplot, color = 'grey', label = 'Inferred cluster')
    ax.bar(np.arange(len(barmeas)), barmeas, color = 'sienna', label = 'Measured cluster')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim([-0.5,len(barplot)-0.5])
    ax.tick_params(labelbottom = False, bottom = False)
    ax.set_ylabel('Number unique clusters')
    if '--symlog' in sys.argv:
        thresh = int(sys.argv[sys.argv.index('--symlog')+1])
        ax.set_yscale('symlog', linthreshy = thresh)
        ax.set_yticks([5,10,20,30,40,50,75,100])
        ax.set_yticklabels([5,10,20,30,40,50,75,100])
    ax.grid( axis = 'y')
    ax.legend()
    
    
    axa = fig.add_subplot(122)
    axa.set_position([0.1,0.05, 0.8,0.43])
    
    shademat = np.zeros(np.shape(clusteroverlaps.T))
    shademat[np.arange(0,len(shademat),2, dtype = int)] = 0.3
    
    if '--upsetplot' in sys.argv:
        axa.imshow(np.ma.masked_array(shademat, mask = clusteroverlaps.T >0), cmap = 'Purples', origin = 'upper', aspect = 'auto', vmin = 0, vmax = 1)
        axa.imshow(np.ma.masked_array(clusteroverlaps.T, mask=clusteroverlaps.T ==0), cmap = 'Greys', origin = 'upper', aspect = 'auto', vmin = 0, vmax = 1)
    
    else:
        print(np.amax(clusterparalogs))
        axa.imshow(np.ma.masked_array(clusterparalogs.T+2, mask=clusteroverlaps.T ==0), cmap = 'Greys', origin = 'upper', aspect = 'auto', vmin = 0, vmax = 5) #np.amax(2+clusterparalogs))
        for cp, cpar in enumerate(clusterparalogs):
            for cq, ca in enumerate(cpar):
                if ca > 0:
                    axa.text(cp, cq, str(np.around(ca,1)), color = 'white', va = 'center', ha = 'center', fontsize = 7)
    
    axa.set_yticks(np.arange(len(clusteroverlaps[0])))
    axa.set_yticklabels(unking)
    axa.set_xlim(ax.get_xlim())
    axa.tick_params(bottom = False)
    axa.spines['top'].set_visible(False)
    axa.spines['left'].set_visible(False)
    axa.spines['bottom'].set_visible(False)
    axa.spines['right'].set_visible(False)
    
    axa.set_xticks(np.arange(len(clusteroverlaps))+0.5, minor = True)
    axa.set_yticks(np.arange(len(clusteroverlaps[0]))+0.5, minor=True)
    axa.grid(which = 'minor', color = 'grey')
    
    axa.set_xticks(np.arange(len(clusteroverlaps)))
    ageclust = []
    for cla in clustage:
        ageclust.append(str(int(cla))+' MYa')
    axa.set_xticklabels(ageclust, rotation = 90)
    
    if '--savefig' in sys.argv:
        print(outname)
        
        fig.savefig(outname+'.svg', bbox_inches = 'tight', dpi = 200)
        figs.savefig(outname+'_speciesclusters.svg', bbox_inches = 'tight', dpi = 200)
        figs2.savefig(outname+'_speciesclusterslost.svg', bbox_inches = 'tight', dpi = 200)
        figs3.savefig(outname+'_speciesevolutionscenario.svg', bbox_inches = 'tight', dpi = 200)
        figs4.savefig(outname+'_specieshomologsequence.svg', bbox_inches = 'tight', dpi = 200)
        figmatev.savefig(outname+'_seqagetospecificityage.svg', bbox_inches = 'tight', dpi = 200)
    plt.show()




