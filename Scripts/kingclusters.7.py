import numpy as np
import sys, os
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import glob
import matplotlib.cm as cm
from sklearn.cluster import DBSCAN
# split KH, RRM
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy import sparse

def assignclade(cluster):
    clades = [['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm', 'Algae', 'Plant', 'Flowering_plant', 'Protists', 'Plasmodium', 'Plant_pathogen', 'Heterokonts'], ['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm', 'Algae', 'Plant', 'Flowering_plant', 'Protists'], ['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm', 'Algae', 'Plant', 'Flowering_plant'], ['Fungi', 'Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm'], ['Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate', 'Insect', 'Fly', 'Worm'], ['Fish', 'Vertebrate', 'Bird', 'Mammal', 'Primate'], [ 'Vertebrate', 'Bird', 'Mammal', 'Primate'], ['Mammal', 'Primate'], ['Insect', 'Fly', 'Worm'], ['Insect', 'Fly'], ['Algae', 'Plant', 'Flowering_plant'], ['Plant', 'Flowering_plant'], [ 'Plant_pathogen', 'Heterokonts'], ['Plasmodium', 'Plant_pathogen', 'Heterokonts']]
    cladenames = ['LECA', 'LECA-some', 'Fungi+Metazoan+Plant', 'Fungi+Metazoan','Metazoan','Vertebrates','Land-vertebrates','Mammals', 'Protostomia', 'Insects', 'Algae+Plants', 'Plants', 'Stramenopiles', 'Alveolata']
    def bestmatch(klist):
        if len(klist)>1:
            fraction = []
            for cla in clades:
                frac = float(len(np.intersect1d(cla,klist)))/float(len(np.union1d(cla, klist)))
                fraction.append(frac)
            return cladenames[np.argmax(fraction)]
        else:
            return klist[0]
    
    assignedclade = []
    for c, clust in enumerate(cluster):
        #print species[clust > 0]
        if np.sum(clust>0) > 1:
            assclust = evodistance[clust>0][0, :]<=np.amax(evodistance[clust>0][:, clust>0])
            asskings = np.unique(kings[assclust])
            #print species[assclust]
            #print asskings
            assignedclade.append(bestmatch(asskings))
        else:
            assignedclade.append(species[clust>0][0])
        #print assignedclade[-1]
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



# function that returns the dictionary with colors for a branch in a linkage_matrix based on the clusters
def make_link_colors(clusters, ulinkage_matrix, colormap):
    unclust, unnum = np.unique(clusters, return_counts = True)
    # set all clusters that have only one leaf to -1
    clusters[np.isin(clusters,unclust[unnum == 1])] = -1
    for c, cl in enumerate(np.unique(clusters)):
        clusters[clusters == cl] = c
    if 1 not in unnum:
        clusters += 1
    curcluster = clusters.astype(float)/np.amax(clusters)
    ccol = {}
    # give clusters a hex color from the provied colormap
    for c, rbg in enumerate(colormap(curcluster)):
        if clusters[c] == 0:
            ccol[c] = 'k'
        else:
            hexi = '#'
            for r in np.array(rbg)[:3]:
                hexi += '%02x' % int(r*255)
            ccol[c] = hexi
    # use the hex colors of clusters to assign colors to branches
    link_cols = {}
    dflt_col = 'k'
    for i, i12 in enumerate(ulinkage_matrix[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(ulinkage_matrix) else ccol[x] for x in i12)
        link_cols[i+1+len(ulinkage_matrix)] = c1 if c1 == c2 else dflt_col
    # return dictionary
    return link_cols



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
    # rank clusters
    ucluster, uclustersize = np.unique(cluster, return_counts = True)
    ucluster = ucluster[-np.argsort(uclustersize)]
    clusterout = np.copy(cluster)
    for c, cl in enumerate(ucluster):
        clusterout[cluster == cl] = c
    cluster = clusterout
    return cluster, clusterdist

# function that returns set of included leafs for each brachpoint in linkage_matrix
def leafset(linkage_matrix):
    leafsets = [[] for i in range(len(linkage_matrix))]
    zmin = len(linkage_matrix) + 1
    for l, link in enumerate(linkage_matrix):
        if link[0] < zmin:
            leafsets[l].append(int(link[0]))
        else:
            leafsets[l] = np.append(leafsets[l],leafsets[int(link[0])-zmin]).astype(int)
        if link[1] < zmin:
            leafsets[l].append(int(link[1]))
        else:
            leafsets[l] = np.append(leafsets[l],leafsets[int(link[1])-zmin]).astype(int)
    return leafsets

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



def plotdendromat(linkage_matrix, matrix, xlabels, link_cols = None, leafnames = None, matrix2 = None, figname = None, dpi = None, vmin = None, vmax = None, vmin2 = None, vmax2 = None, grid = False, grid2 = False, colormap = None, colormap2 = None, line = None):
    if dpi is None:
        dpi = 100
    if matrix2 is None:
        fig = plt.figure(figsize = (1+0.2*len(matrix[0]),8), dpi = dpi)
        nmat = 2
    else:
        fig = plt.figure(figsize = (1+0.2*len(matrix[0]),12), dpi = dpi)
        nmat = 3
    if link_cols is None:
        link_cols = {}
        for i in range(len(linkage_matrix)):
            link_cols[i+1+len(linkage_matrix)] = 'k'
            
    if leafnames is None:
        leafnames = ['' for i in range(len(matrix[0]))]
        ax1.tick_params(bottom = False, labelbottom = False)
        
    ax1 = fig.add_subplot(nmat, 1,1)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    dn = dendrogram(linkage_matrix, ax = ax1, link_color_func = lambda k: link_cols[k], labels = leafnames, leaf_rotation = 90, leaf_font_size = 6)
    
    if line is not None:
        lim = ax1.get_xlim()
        ax1.plot(lim, [line,line], 'red', ls = '--')
    
    sort = dn['leaves']
    ax2 = fig.add_subplot(nmat, 1,2)
    if len(matrix) == len(sort):
        matrix = matrix[sort][:,sort]
    else:
        matrix = matrix[:,sort]
    
    if colormap is None:
        colormap = 'Greys'
    ax2.imshow(matrix, aspect = 'auto', cmap = colormap, vmin = vmin, vmax = vmax)
    ax2.tick_params(left = False, labelleft = False, top = True)
    if matrix2 is not None:
        ax2.tick_params(bottom = False, labelbottom = False)
    else:
        ax2.set_xticks(np.arange(len(matrix[0])))
        ax2.set_xticklabels(xlabels[sort], rotation = 90)
    if grid2:
        ax2.set_xticks(np.arange(len(matrix[0]))-0.5, minor = True)
        ax2.set_yticks(np.arange(len(matrix))-0.5, minor = True)
        ax2.grid(color ='grey', which = 'minor')
    
    if matrix2 is not None:
        ax3 = fig.add_subplot(nmat, 1,3)
        if len(matrix2) == len(sort):
            matrix2 = matrix2[sort][:,sort]
        else:
            matrix2 = matrix2[:,sort]
        if colormap2 is None:
            colormap2 = 'Greys_r'
        ax3.imshow(matrix2, aspect = 'auto', cmap = colormap2, vmin = vmin, vmax = vmax)
        ax3.tick_params(left = False, labelleft = False)
        ax3.set_xticks(np.arange(len(matrix2[0])))
        ax3.set_xticklabels(xlabels[sort], rotation = 90)
        if grid2:
            ax3.set_xticks(np.arange(len(matrix2[0]))-0.5, minor = True)
            ax3.set_yticks(np.arange(len(matrix2))-0.5, minor = True)
            ax3.grid(color ='grey', which = 'minor')
    
    if figname is not None:
        fig.savefig(figname, bbox_inches = 'tight')
        print 'SAVED', figname
    else:
        fig.tight_layout()
        plt.show()
    plt.close()
    return


def printcluster(cluster, names, target):
    cls = np.unique(cluster)
    for c, cl in enumerate(cls):
        if target in names[cluster == cl]:
            print names[cluster == cl]



# find cores by triangles and min flatdist
# refine cores with agglomerative clustering
# extend cores to orthologs with < confcut
# extetnd new set to homologs if species tree is not extended

# function uses idmat to cluster with scipy linkage, using method as linkage
def agglomerative2(hidmat, flatdist, confcut, correctmat = None, maskmat = None, method = 'average', plot_all = False, **kwargs):
    
    # combine all RBRs that are triplet connected and similar binding according to flatdist
    cluster = clusterall(hidmat*(flatdist<=confcut))
    cls, cln = np.unique(cluster, return_counts = True)
    print 'Clusterall', len(np.unique(cluster)), np.unique(cln, return_counts = True)
    
    # resort clusters by size
    ncluster = -np.ones(len(cluster))
    cls = cls[np.argsort(-cln)]
    for c, cl in enumerate(cls):
        ncluster[cluster == cl] = c
    cluster = ncluster
    maxcl = np.amax(cluster)
   
    printcluster(cluster, realprotnames, 'Syp__RRM-RRM-RRM')
    printcluster(cluster, species[protspecies], 'A1CF')
    
    for c in range(len(cls)):
        cl = cls[c]
        cmask = np.where(cluster == cl)[0]
        if len(cmask) > 1:
            cmat = flatdist[cmask][:,cmask]
            if (cmat > confcut).any():
                linkage_matrix = linkage(cmat[np.triu_indices(len(cmat),1)], method=method)
                if (linkage_matrix[:,2]>confcut).any():
                    if correctmat is not None:
                        correctmask = (correctmat[cmask][:,cmask])&(cmat< confcut)
                        cmat[correctmask] /= 2.
                        linkage_matrix = linkage(cmat[np.triu_indices(len(cmat),1)], method=method)
                    nc, ncdist = splittree(linkage_matrix, confcut)
                    for n in np.unique(nc)[1:]:
                        maxcl +=1
                        cluster[cmask[nc == n]] = maxcl
                    if len(cmask) < 500 and plot_all:
                        if 'figname' in kwargs:
                            figname = kwargs['figname']+'_precluster_'+method+str(c)+'.jpg'
                        else:
                            figname = outname+'_precluster_'+method+str(c)+'.jpg'
                        plotdendromat(linkage_matrix, cmat, species[protspecies[cmask]], link_cols = make_link_colors(cluster[cmask], linkage_matrix, cm.gist_rainbow), leafnames = cluster[cmask].astype(int), matrix2 = None, figname = figname, dpi = None, vmin = 0, vmax = 0.4, vmin2 = None, vmax2 = None, colormap = 'Greys_r', colormap2 = None, grid = False, grid2 = False)
    
    cls, cln = np.unique(cluster, return_counts = True)
    print len(np.unique(cluster)), np.unique(cln, return_counts = True)
    
    ncluster = -np.ones(len(cluster))
    cls = cls[np.argsort(-cln)]
    for c, cl in enumerate(cls):
        ncluster[cluster == cl] = c
    cluster = ncluster.astype(int)
    maxcl = np.amax(cluster)
    
    printcluster(cluster, realprotnames, 'Syp__RRM-RRM-RRM')
    printcluster(cluster, species[protspecies], 'A1CF')
    
    if maskmat is not None:
        cls, cln = np.unique(cluster, return_counts = True)
        # unite clusters within species if their assigned clusters share mean less than cutoff
        for c, cl in enumerate(cls):
            smask = np.where(cluster == cl)[0]
            potspecies = np.unique(protspecies[smask])
            if len(potspecies) > 0:
                potprot = []
                for ps in potspecies:
                    spmask = smask[protspecies[smask] == ps]
                    canhavelink = flatdist[spmask]
                    canhavelink[~maskmat[spmask]] = confcut+ 0.1
                    havelink = np.where((np.sum(canhavelink <= confcut,axis = 0)>0)&(protspecies == ps))[0]
                    potprot.append(havelink)
                potclust = cluster[np.concatenate(potprot)]
                potclust = potclust[potclust!= cl]
                potclust = np.unique(potclust)
                if len(potclust) > 0:
                    potprot = np.where(np.isin(cluster, potclust))[0]
                    cmask = np.union1d(smask, potprot)
                    cmat = flatdist[cmask][:,cmask]
                    if correctmat is not None:
                        correctmask = (correctmat[cmask][:,cmask])&(cmat< confcut)
                        cmat[correctmask] /= 2.
                    linkage_matrix = linkage(cmat[np.triu_indices(len(cmat),1)], method=method)
                    nc, ncdist = splittree(linkage_matrix, confcut)
                    ncwith = nc[cluster[cmask] == cl]
                    nctotest, nctestcount = np.unique(cluster[cmask[np.isin(nc,ncwith)]], return_counts = True)
                    nctotest = nctotest[-np.argsort(nctestcount)]
                    nctotest = nctotest[nctotest != cl]
                    if len(nctotest)>0:
                        if len(cmask) < 500 and plot_all:
                            
                            if 'figname' in kwargs:
                                figname = kwargs['figname']+'_postcluster_'+method+str(c)+'.jpg'
                            else:
                                figname = outname+'_postcluster_'+method+str(c)+'.jpg'
                            showmat = hidmat[cmask][:,cmask].astype(int) + ((flatdist[cmask][:,cmask] <= confcut)*(hidmat[cmask][:,cmask])).astype(int)
                            plotdendromat(linkage_matrix, showmat, species[protspecies[cmask]], link_cols = make_link_colors(nc, linkage_matrix, cm.gist_rainbow), leafnames = cluster[cmask].astype(int), matrix2 = None, figname = figname, dpi = None, vmin = 0, vmax = 2, vmin2 = None, vmax2 = None, colormap = 'Greys', colormap2 = None, grid = False, grid2 = False)
                        cluster[np.isin(cluster,nctotest)] = cl
                        
    cls, cln = np.unique(cluster, return_counts = True)
    print len(np.unique(cluster)), np.unique(cln, return_counts = True)
    
    ncluster = -np.ones(len(cluster))
    cls = cls[np.argsort(-cln)]
    for c, cl in enumerate(cls):
        ncluster[cluster == cl] = c
    cluster = ncluster.astype(int)
    maxcl = np.amax(cluster)

    printcluster(cluster, realprotnames, 'Syp__RRM-RRM-RRM')
    printcluster(cluster, species[protspecies], 'A1CF')
    

    if plot_all:
        # look at subclusters that include several clusters
        if correctmat is not None:
            flatdist[correctmat & (flatdist<confcut)] /= 2. 
        linkage_matrix = linkage(flatdist[np.triu_indices(len(flatdist),1)], method=method)
        fcluster, fclusterdist = splittree(linkage_matrix, cutlink = 2.5*confcut)
        ufcluster, ufclnum = np.unique(fcluster, return_counts = True)

        for ufc in ufcluster:
            ufmask = np.isin(cluster, cluster[fcluster == ufc])
            ufmat = flatdist[ufmask][:, ufmask]
            if 'figname' in kwargs:
                figname = kwargs['figname']+'_allcluster_'+method+str(ufc)+'.jpg'
            else:
                figname = outname+'_allcluster_'+method+str(ufc)+'.jpg'
            if len(np.unique(cluster[ufmask]))>2 and len(ufmat) < 500:
                linkage_matrix = linkage(ufmat[np.triu_indices(len(ufmat),1)], method=method)
                showmat =  (flatdist[ufmask][:,ufmask] <= confcut).astype(int) + 0.5* hidmat[ufmask][:,ufmask].astype(int) 
                link_cols = make_link_colors(cluster[ufmask], linkage_matrix, cm.gist_rainbow)
                plotdendromat(linkage_matrix, showmat, species[protspecies[ufmask]], link_cols =link_cols, leafnames = cluster[ufmask].astype(int), matrix2 = None, figname = figname, dpi = None, vmin = 0, vmax = 2, vmin2 = None, vmax2 = None, colormap = 'Greys', colormap2 = None, grid = False, grid2 = False, line = confcut)
                print 'SAVED', figname

    return cluster
    
# false negatives are defined as clusters with a single element but have more than uniteval identity to another cluster
def cleanfncluster(idmat, clusterid, uniteval, cluster, maskmat = None):
    uniquecluster, uniqueclen = np.unique(clusterid, return_counts = True)
    idmat = np.copy(idmat)
    np.fill_diagonal(idmat, 0)
    for c, cl in enumerate(uniquecluster):
        # check if cl is single cluster
        if (uniqueclen[c] == 1):
            # check for largest similar data point
            if maskmat is not None:
                unitein = np.where(maskmat[clusterid == cl][-1])[0]
                unitein = unitein[np.argmax(idmat[clusterid == cl][-1][unitein])]
            else:
                unitein = np.argmax(idmat[clusterid == cl][-1])
            # check if best data point fullfills uniteval condition
            if idmat[unitein, clusterid == cl][0] >= uniteval:
                # check if 
                core = np.amax(evodistance[cluster[clusterid[unitein]]>0][:, cluster[clusterid[unitein]]>0])
                evosort = np.unique(evodistance[cluster[clusterid[unitein]]>0][0])
                evosort = evosort[evosort>core]
                if len(evosort) > 0:
                    core = evosort[0]
                core = evodistance[cluster[clusterid[unitein]]>0][0] <= core
                if np.sum(core*cluster[c]) != 0:
                    #print 
                    #print ''.join((cluster[cl]>0).astype(int).astype(str))
                    #print ''.join((cluster[clusterid[unitein]]>0).astype(int).astype(str))
                    #print ''.join(core.astype(int).astype(str))
                    clusterid[clusterid == cl] = clusterid[unitein]
    # reset cluster identifiers
    clusteridout = np.copy(clusterid)
    for c, cl in enumerate(np.unique(clusterid)):
        clusteridout[clusterid == cl] = c
    clusterid = clusteridout
    
    clusterassign, unnum = np.unique(clusterid, return_counts = True)
    cluster = []
    clusterage = []

    for c, cl in enumerate(clusterassign):
        clspecies = protspecies[clusterid == cl]
        clspecies, clspecnum = np.unique(clspecies, return_counts = True)
        clpat = np.zeros(len(species), dtype = np.int8)
        clpat[clspecies] = clspecnum
        cluster.append(clpat)
        clusterage.append(np.amax(evodistance[clspecies][:, clspecies]))
    cluster = np.array(cluster, dtype = np.int8)
    clusterage = np.array(clusterage)
    
    return cluster, clusterage, clusterid



# False positives can be spotted as sequences that are not one-to-one orthologs to a member of the cluster and would result in younger common ancestor if removed
def cleanfpcluster(evidmat, evodistance, clusterage, cluster, clusterid, pspecies):
    cmax = np.amax(clusterid)
    print 'FP clean before', len(cluster)
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
                            #print c
                            #print ''.join(cspec.astype(int).astype(str))
                            cmax += 1
                            clusterid[cmask[nto]] = cmax
                            cluster = np.append(cluster,[np.isin(arange,espec[nto]).astype(int)] ,axis = 0)
                            clusterage = np.append(clusterage, [0])
                            cluster[c, espec[nto]] = 0
                            cspec = cluster[c] > 0
                            clusterage[c] = ncluage
                            repeat = False
                            break
        
        #if np.sum(cspec.astype(float))/np.sum(coverage.astype(float)) <0.1:
            
            #print ''.join((cluster[c]>0).astype(int).astype(str))
            #print species[protspecies[cmask]], notonetoone
            #print evidmat[cmask][:,cmask]
        
        if repeat:
            c += 1
            
        if c == len(cluster) - 1:
            break
    print 'FP clean after', len(cluster)
    return cluster, clusterage, clusterid
                       

# compute homolog matrix, 1-1 ortholog matrix and triangle 1-1 ortholog matrix pairs
# these matrices contain additional information about the sequence identity relationship between species
def homologmatrix(idmat, pspecies, testortholog = False, testtriangle = False, maskmat = None, include_self = False):
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


# for each cluster, use sequence identity or homolog matrix to identify most likely common ancestor of cluster
# and determine their relationship
def determine_family(clusters, clusterage, evodistance, clusterid, simmat, trust = None, maskmat = None, maxin = 1, maxo = 0):
    familyage = np.zeros(len(clusters))
    family = -np.ones(len(clusters))
    famevidence = np.zeros(len(clusters))
    family_relation = np.zeros(len(clusters), dtype = np.int8)
    
    for p in range(len(clusters)):
        if p % 100 == 0:
            print p
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
            #print paralogage,  sortevo[0], species[insidespec], species[outsidespec]
            fstages = []
            for mi in range(maxout+1):
                considerstages = evodistance[outsidespec]
                for co, cost in enumerate(considerstages):
                    if sortevo[mi] in cost:
                        cost[cost > sortevo[mi]] = sortevo[mi]
                        #print mi, species[outsidespec[co]], sortevo[mi], np.unique(cost), np.unique(cost)[-1-maxin:]
                        fstages.append(np.unique(cost)[-1-maxin:])
            fstages = np.unique(np.concatenate(fstages))
            fstages = fstages[fstages != 0]
            fstages = np.unique(np.append(fstages, sortevo))
            #print fstages
            #print np.unique(clusterage[potential_family])
            potential_family = potential_family[np.isin(clusterage[potential_family], fstages)]
            
            #print np.unique(clusterage[potential_family])
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
            
            #print np.unique(clusterage[potential_family])
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
        
        
     
    


if __name__ == '__main__':
    
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
        measprots = latentmeas['names']
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
        protsimnames = protein_simfile['names']
        
        # maximal number of domains that are considered to be in one protein sequence
        maxdom = 5
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
                protdomains[lok,p] = 1
                protdomains2[-1+min(maxdom,dtnum[list(dtypes).index('RRM')]), p] = 1
                
            if 'KH' in dtypes:
                lok = [maxdom+min(max(0,x),maxdom-1) for x in range(-2+min(maxdom+1,dtnum[list(dtypes).index('KH')]), 1+min(maxdom+1,dtnum[list(dtypes).index('KH')]))]
                protdomains[lok,p] = 1
                protdomains2[maxdom-1+min(maxdom,dtnum[list(dtypes).index('KH')]), p] = 1
            
            domaincounts[p] = np.sum(dtnum)
            realprotnames[p] = psname[-2]+'__'+psname[-1]
        
        # remove human pseudo proteins
        protsim = protsim[remove][:, remove]
        protsimnames = protsimnames[remove]
        realprotnames = realprotnames[remove]
        protdomains = protdomains[:, remove]
        protdomains2 = protdomains2[:, remove]
        domaincounts = domaincounts[remove]
        
        if '--correct_domaincount' in sys.argv:
            print 'Correcting SeqID with domain count'
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
            latentnames.append(slat['names'])
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
            
            print 'Correct for minimum sequence length', lcorlength
            print len(latentnames)
            latentnames = latentnames[sort]
            protspecies = protspecies[sort]
            latent = latent[sort]
            protsimnames = protsimnames[sort]
            realprotnames = realprotnames[sort]
            protsim = protsim[sort][:, sort]
            protdomains = protdomains[:, sort]
            protdomains2 = protdomains2[:, sort]
            print len(latentnames)

        # generate 3 matrices that mask only rbps that can be compared depending on their domain composition
        # has ones for two RBPs that have equal number of domains or one difference
        protdommat = np.dot(protdomains2.T, protdomains) > 0
        # has ones only if RBPs have same number of RBDs
        protdomains = np.array([np.sum(protdomains[:maxdom], axis = 0) > 0, np.sum(protdomains[maxdom:], axis = 0) > 0])
        protdommat2 = np.dot(protdomains.T, protdomains) > 0
        # can have ones if x1-1 = x2+1
        protdommat3 = np.dot(protdomains2.T, protdomains2) > 0
        
        # set identities of rbps that cannot be compared to 0
        #protsim[~protdommat2] = 0.
        
        # compute latent distances to measured rbps
        latdist = cdist(measlatent, latent, 'cosine')
        closestdist = np.amin(latdist, axis = 0)
        closestmeas = measprots[np.argmin(latdist, axis = 0)]
        # compute latent distances between rbps
        flatdist = cdist(latent, latent, 'cosine')
        
        # compute homolog and ortholog matrix
        hmatrix, omatrix, tmatrix = homologmatrix(protsim, protspecies, testortholog = True, testtriangle= True, maskmat = protdommat)
        
        # correctneighbor assigns lower cosine values to RBPs if their ortholog in the most related species is closest to each other in latent space.
        if '--correctneighbor' in sys.argv:
            evoneigbor = evodistance[protspecies]
            evobest = np.sort(evoneigbor, axis = 1)[:,1]
            evoneigbor = evoneigbor[:,protspecies]
            evoneigbor = evoneigbor-evobest == 0
            evoneigbor = evoneigbor.astype(int)
            evoneigbor = (evoneigbor + evoneigbor.T) > 0 
            evoneigbor = omatrix * evoneigbor
        else:
            evoneigbor = None
        
        
        #perform agglomerative clustering with linkage
        # use protdommat which indicates if rbps with same domaintype can be assinged into same cluster
        clusterid = agglomerative2(hmatrix, flatdist, confcut, correctmat = evoneigbor, maskmat = protdommat)
        clusterassign, unnum = np.unique(clusterid, return_counts = True)
        print 'Unique cluster detected', len(clusterassign), np.amax(clusterid)
        print 'Size distribution', np.unique(unnum, return_counts = True)
        
        cluster = []
        clusterage = []
        for c, cl in enumerate(clusterassign):
            clspecies = protspecies[clusterid == cl]
            clspecies, clspecnum = np.unique(clspecies, return_counts = True)
            clpat = np.zeros(len(species), dtype = np.int8)
            clpat[clspecies] = clspecnum
            cluster.append(clpat)
            clusterage.append(np.amax(evodistance[clspecies][:, clspecies]))
        cluster = np.array(cluster, dtype = np.int8)
        clusterage = np.array(clusterage)
        
        # remove false positives from cluster if they harbour species outside the core of the cluster and have no support of being a one-to-one ortholog
        cluster, clusterage, clusterid = cleanfpcluster(tmatrix, evodistance, clusterage, cluster, clusterid, protspecies)        
        print 'Unique cluster after FP correction', len(cluster), np.amax(clusterid), len(np.unique(clusterid))
        
        # care for false negatives, which are all rbps that are 70% identical and also have same domain structure. 
        cluster, clusterage, clusterid  = cleanfncluster(protsim, clusterid, psimcut, cluster, maskmat = protdommat3)

        
        
        print np.unique(clusterage, return_counts = True)
        clusterassign, unnum = np.unique(clusterid, return_counts = True)
        print 'Unique cluster after FN correction', len(clusterassign), np.unique(unnum, return_counts = True)
        
        # determine family relationship choosing from different methods
        if '--allidlink' in sys.argv:
            outname += '_idlink'
            family_relation, familyage, family, famevidence = determine_family(cluster, clusterage, evodistance, clusterid, protsim, trust= plowcut)
        elif '--homolink' in sys.argv:
            outname += '_homolink'
            family_relation, familyage, family, famevidence = determine_family(cluster, clusterage, evodistance, clusterid, (protsim+protsim*(tmatrix.astype(float)+omatrix.astype(float))/2.)/2., trust= None)
        else:
            family_relation, familyage, family, famevidence = determine_family(cluster, clusterage, evodistance, clusterid, protsim, trust= plowcut, maskmat = hmatrix)
        
        print 'Final age distribution of RSSGs', np.unique(family_relation, return_counts = True)
        
        
        clustermeasbool = np.isin(clusterid, clusterid[closestdist < confcut])
        
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
        
        print 'SAVED', outname+'.txt'
        np.savetxt(outname+'.txt', np.append(np.array([latentnames, realprotnames, clusterid, clustermeasbool.astype(int), np.around(closestdist,3), closestmeas, np.sum(clusterprot>=1,axis = 1).astype(int), clusterages, seqage, evolutionorigin, evolutionids, evolutiontype, protspecies]).T, clusterprot, axis = 1), header = 'ProtId Protname Cluster Measured Closest_jple CloestRBR N_species Age_spec Age_seq Parent ParentSID Evolution SpeciesId; '+' '.join(species), fmt = '%s')
        
        # SAVE text file that is only for clusters and give one to three RBP names for each cluster: Names from certain species are preferred.
        speciesrank = ['Homo_sapiens', 'Mus_musculus', 'Macaca_mulatta', 'Gallus_gallus', 'Danio_rerio', 'Xenopus_tropicalis', 'Saccharomyces_cerevisiae', 'Schizosaccharomyces_pombe', 'Drosophila_melanogaster', 'Drosophila_ananassae', 'Caenorhabditis_elegans', 'Caenorhabditis_briggsae', 'Arabidopsis_thaliana', 'Arabidopsis_lyrata', 'Zea_mays', 'Tetraodon_nigroviridis', 'Oryza_sativa', 'Canis_familiaris', 'Latimeria_chalumnae', 'Geospiza_fortis', 'Musca_domestica']
        nameranking = np.zeros(len(protspecies))
        for s, spec in enumerate(species):
            if spec in speciesrank:
                nameranking[protspecies == s] = speciesrank.index(spec)
            else:
                nameranking[protspecies == s] = len(species)
            
        
        # find cluster representatives
        # and best closest measurement for cluster
        repclustnames = []
        bestclustmeasurement = []
        repdomains = []
        for c, cl in enumerate(clusterassign):
            exnames = realprotnames[clusterid == c][np.argsort(nameranking[clusterid == c])]
            exdoms = np.copy(exnames)
            for e, en in enumerate(exnames):
                exnames[e], exdoms[e] = en.rsplit('__',1)
            exnames = exnames[:3]
            exdoms, edomnum = np.unique(exdoms,return_counts = True)
            exdoms = exdoms[np.argmax(edomnum)]
            exmeas, exmeasn = np.unique(closestmeas[clusterid == c], return_counts = True)
            exmeas = exmeas[np.argmax(exmeasn)]
            repclustnames.append(','.join(exnames))
            bestclustmeasurement.append(exmeas)
            repdomains.append(exdoms)
        

        
        assignedclade = assignclade(cluster)
        berelative, berelativenum = np.unique(family.astype(int), return_counts = True)
        isrelative = np.zeros(len(family), dtype = int)
        isrelative[berelative] = berelativenum
        
        clustermeasbool = np.isin(np.arange(len(cluster)),np.unique(clusterid[clustermeasbool]))
        dataarray = np.array([clusterassign, clustermeasbool, np.array(repclustnames), np.array(repdomains), np.array(bestclustmeasurement), np.sum(cluster,axis = 1).astype(int), clusterage, familyage, family.astype(int), np.around(famevidence,0), family_relation, np.sum(cluster>0,axis = 1), np.array(assignedclade), isrelative]).T
        np.savetxt(outname+'_cluster.txt', np.append(dataarray, cluster, axis = 1), header = 'ClusterId;Measured;Proteinnames;Domainclass;RNAcmpt;Num_RBPs;Age_spec;Age_seq;RelativeRSSG;RelativeSID;Evolution;Num_species;Clade;Num_ancestor;'+';'.join(species), fmt = '%s', delimiter = ';')
        
        uniquecl, uniqueclindex = np.unique(clusterid, return_index = True)
        seqage = seqage[uniqueclindex] 
        clustermeasured = cluster[clustermeasbool]
        evolutiontail = family_relation
        sys.exit()
        
    elif '--datasheet' in sys.argv:
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
            outname = sys.argv[sys.argv.index('--outname')
        else:
            outname = os.path.splitext(sys.argv[sys.argv.index('--datasheet')+1])[0]
        cluster = clusterprot[:,12:].astype(int)[clusterindex]
        closestdist = clusterprot[:,3].astype(float)
        closestmeasured = clusterprot[:,4]
        clustermeasured = cluster[clustermeasbool]
        print 'Measured', len(clustermeasured)
        clusterlen = clusterprot[:,5].astype(int)
        clusterage = clusterprot[:,6].astype(float)[clusterindex]
        seqage = clusterprot[:,7].astype(float)[clusterindex]
        protspecies = clusterprot[:,11].astype(int)
        evolutiontype = clusterprot[:,10].astype(int)
        evolutiontail = evolutiontype[clusterindex]
        
    else:
        print 'Provide data sheet or analyze latent space'
        sys.exit()
    
    
    
    if '--confident' in sys.argv:
        confident = float(sys.argv[sys.argv.index('--confident')+1])
        
        confidentclusters = np.unique(clusterid[closestdist <= confident])
        nonconfidentclusters = clusterassign[~np.isin(clusterassign, confidentclusters)]
        
        cid, cidd = np.unique(clusterid, return_index = True)
        clustlen = clusterlen[cidd]
        print 'Confident Noncofident'
        print 'Mean', np.mean(clustlen[confidentclusters]), np.mean(clustlen[nonconfidentclusters])
        print 'Number', len(confidentclusters), len(nonconfidentclusters)
        
        if '--nonconfident' in sys.argv:
            outname += '_nonconfidence'+str(confident)
            nonconfidentclusters, confidentclusters = confidentclusters, nonconfidentclusters
        else:
            outname += '_confidence'+str(confident)
        

    
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
        print domainspec
        print 'Mean', np.mean(clustlen[confidentclusters]), np.mean(clustlen[nonconfidentclusters])
        print 'Number', len(confidentclusters), len(nonconfidentclusters)
        outname += '_'+domainspec
    
    
    if '--domainspecific' in sys.argv or '--confident' in sys.argv:
        cluster = cluster[np.isin(clusterassign,confidentclusters)]
        clustermeasbool = clustermeasbool[np.isin(clusterassign,confidentclusters)]
        clustercenter = clustercenter[np.isin(clusterid,confidentclusters)]
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
    
    unking = []
    for kun in kings:
        if kun not in unking:
            unking.append(kun)
    unking = np.array(unking)
    
    
    uniqueage, ageid = np.unique(clusterage, return_index = True)
    uniqueage, ageid = uniqueage[::-1], ageid[::-1]
        
    shortspec = []
    for s, spec in enumerate(species):
        shortspec.append(spec[0])
    shortspec = np.array(shortspec)
    

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
            #print '\n', uniage
            #print ''.join(cladespec.astype(int).astype(str))
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
                #print compage
                for gi, ga in enumerate(gains):
                    #print ''.join(cladespec.astype(int).astype(str)), ''.join(gaininclude[gi].astype(int).astype(str))
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
    
    # generate subclusters based on distribution of age clusters
    kingcluster = np.zeros((len(cluster), len(unking)))
    for spco in np.unique(speciscolors)[:-1]:
        spcowhere = np.where(speciscolors == spco)[0]
        for c, csp in enumerate(cluster[spcowhere]):
            ckings = np.unique(kings[csp>=1])
            for ck in ckings:
                kingcount = csp[kings == ck]
                kingcluster[spcowhere[c], list(unking).index(ck)] = np.mean(kingcount[kingcount>0])
            #print kingcluster[spcowhere[c]], csp
        
    kingclustermeasured = kingcluster[clustermeasbool]
    
    

    if '--plot_branchpoints' in sys.argv:
        total = []
        for x, xname in enumerate(xnames):
            total.append(np.sum(lossspec[x]))
            total.append(np.sum(gainsplit[x]))
            #print maintains[x], np.sum(maintains[x])
        #print np.amax(total), total
        for x, xname in enumerate(xnames):
            #print x, xname, gains[x], losses[x], maintains[x]
            fig =plt.figure(figsize = (4.5,2.5), dpi = 250)
            ax = fig.add_subplot(211)
            ax.set_position([0.35, 0.5, 0.3, 0.3])
            #ax.set_ylim([0, np.amax(total)])
            ax.set_title(xname)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_linewidth(1.4)
            #ax.spines['right'].set_visible(False)
            ax.tick_params(labelright = False, labelleft = False, left = False, labelbottom = False, bottom = False, width = 1.4)
            #ax.barh([3.], np.cumsum(maintains[x])[::-1], color = colors[:x+1][::-1], lw = 0. , edgecolor = 'grey')
            ax.barh([2.5], np.cumsum(evolution[x])[::-1], color = cm.tab20([15, 8, 6, 10, 12]), lw = 0. , edgecolor = 'grey')
            ax.barh([1.5], np.cumsum(gainsplit[x])[::-1], color = colors[:x+1][::-1], lw = 0. , edgecolor = 'grey')
            lossx = np.cumsum(lossspec[x])[::-1]
            ax.barh([0], lossx, color = colors[:x+1][::-1])
            #ax.set_yticks([1])
            #ax.text(np.sum(evolution[x])+1,3, str(int(np.sum(evolution[x]))), ha = 'left', va = 'center')
            ax.text(np.amax(lossx)+1,0, str(-int(np.amax(lossx))), ha = 'left', va = 'center')
            ax.text(np.sum(gainsplit[x])+1,1.5, str(int(np.sum(gainsplit[x]))), ha = 'left', va = 'center')
            ax.set_ylim([-0.8, 3.3])
            ax.set_xlim([0, np.amax(total)])
            #ax.set_yticklabels(['Gains', 'Losses'], rotation = 60, ha = 'right', va = 'top')
            
            axp = fig.add_subplot(121)
            axp.set_position([0.2, 0.55, 0.12*2.5/4.5, 0.12])
            axp.set_axis_off()
            #circle = plt.Circle((0.,0.), np.sqrt(gains[x]/np.pi), facecolor = colors[x], ls = '-', lw = 0., edgecolor = 'k')
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
            print outname+xname.replace(' ', "_")+'.png'
            plt.close()
        
    
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
        # compuate all clusterages that are in the 0 color categorie
        addcolor, addcolornum = np.unique(-clusterage[(cluster[:,s] >=1)*(speciscolors == uncolor[-1])], return_counts = True)
        # collect colors for lost specificities, don't have the color of single ones
        species_colorslost.append(np.array(colors)[uncolor[:-1].astype(int)][::-1])
        # add a second color for overlap with nearest neighbor species
        morenum = np.append(unnum[:-1], addcolornum)
        morecolor = np.array(colors)[uncolor.astype(int)][:-1]
        for a, adc in enumerate(addcolor):
            changecolor = np.array(colors[-1])
            changecolor[-1] = 0.8 + 0.2*-adc/max(-addcolor[0],1.)
            morecolor = np.append(morecolor, [changecolor], axis = 0)
        
        speciesbars.append(np.cumsum(morenum)[::-1])
        species_colors.append(morecolor[::-1])
        species_age.append(clusterage[cluster[:,s] >= 1][unind])
        
        #print speciesbars[-1]
        
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
        #print spec, potages
        specorn = np.zeros((len(potages), len(potages)))
        seqtocluster, staing, seqtoclustnum = np.unique(zip(seqage[cluster[:,s]>=1], clusterage[cluster[:,s]>=1]), return_counts = True, return_index= True,  axis = 0)
        for tc, stc in enumerate(seqtocluster): 
            #print stc, seqtoclustnum[tc]
            specorn[potages.index(stc[0]), potages.index(stc[1])] = seqtoclustnum[tc]
        speciesorigin.append(specorn)
        
        unnumspecies = totalnum[np.isin(totalcluster, uncolor)]-unnum 
        speciesbarslost.append(np.cumsum(unnumspecies[:-1])[::-1])
        sevot, sevotn = np.unique(evolutiontail[cluster[:, s]>=1], return_counts = True)
        evoln = np.zeros(5)
        evoln[sevot] = sevotn
        speciesevolution.append(np.cumsum(evoln)[::-1])
        #print speciesevolution[-1]
 
    
    
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
    
    if '--horizontalbar' in sys.argv:
        figs2 = plt.figure(figsize = (len(species)*0.28, 4), dpi = 100)
    else:
        figs2 = plt.figure(figsize = (4.5, len(species)*0.28), dpi = 100)
    axs2 = figs2.add_subplot(111)
    for s, spec in enumerate(species):
        if '--horizontalbar' in sys.argv:
            axs2.bar(np.ones(len(speciesbarslost[s]))*s, speciesbarslost[s], color = species_colorslost[s], linewidth = 0., edgecolor = 'k')
        else:
            #print speciesbarslost[s]
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
    
    # percentage
    if '--horizontalbar' in sys.argv:
        figs3 = plt.figure(figsize = (len(species)*0.28, 4), dpi = 100)
    else:
        if '--evopercentage' in sys.argv:
            figs3 = plt.figure(figsize = (3., len(species)*0.28), dpi = 100)
        else:
            figs3 = plt.figure(figsize = (4.5, len(species)*0.28), dpi = 100)
            
    axs3 = figs3.add_subplot(111)
    for s, spec in enumerate(species):
        if '--evopercentage' in sys.argv:
            #print speciesevolution[s]
            speciesevolution[s] = speciesevolution[s] - speciesevolution[s][-1]
            speciesevolution[s][speciesevolution[s]<0] = 0
            speciesevolution[s] = speciesevolution[s]*100./max(1,np.amax(speciesevolution[s]))
            
        if '--horizontalbar' in sys.argv:
            axs3.bar(np.ones(len(speciesevolution[s]))*s, speciesevolution[s], color = cm.tab20([15, 8,6,10,12]), linewidth = 0., edgecolor = 'k')
        else:
            axs3.barh(np.ones(len(speciesevolution[s]))*-s, speciesevolution[s], color =cm.tab20([15, 8,6,10,12]), linewidth = 0., edgecolor = 'k')
    axs3.spines['top'].set_visible(False)
    axs3.spines['right'].set_visible(False)
    if '--horizontalbar' in sys.argv:
        axs3.set_xticks(np.arange(len(species)))
        axs3.set_xticklabels(species, rotation = 90)
    else:
        axs3.set_yticks(-np.arange(len(species)))
        axs3.set_yticklabels(species, rotation = 0)
        axs3.set_ylim([-len(species), .5])
        print np.amax(np.concatenate(speciesevolution))
        axs3.set_xticks(np.arange(0, np.amax(np.concatenate(speciesevolution)), 20, dtype = int))
        axs3.grid(axis = 'x')
    
    if '--evopercentage' in sys.argv:
        axs3.set_xlabel('Percentage evolutionary relationship')
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
        

    
    cladecolor = ['firebrick', 'blue', 'skyblue', 'k', 'red', 'sienna', 'green', 'limegreen', 'orange', 'grey', 'seagreen', 'purple', 'goldenrod']
    chosenspecies = ['Schizosaccharomyces_pombe', 'Danio_rerio', 'Gallus_gallus', 'Homo_sapiens', 'Drosophila_melanogaster', 'Caenorhabditis_elegans', 'Selaginella_moellendorffii', 'Arabidopsis_thaliana', 'Trichomonas_vaginalis', 'Plasmodium_falciparum', 'Albugo_laibachii',  'Thalassiosira_pseudonana']
    if '--separate_timeplot' not in sys.argv:
        fige = plt.figure(figsize = (5, 5), dpi = 200)
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
            axe.plot(-species_age[s][:-1], -speciesbarslost[s][:-1], color = cladecolor[q], marker = 'v', markersize = 3., ls= '--', label = spec)
        
        fige.savefig(outname+'_species_clusterevolution.jpg', bbox_inches = 'tight', dpi = 300)
        

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
        
        fige.savefig(outname+'_specificitytiming.jpg', bbox_inches = 'tight', dpi = 300)
        
        
    else:
        for q, spec in enumerate(chosenspecies):
            fige = plt.figure(figsize = (6, 3), dpi = 200)
            axe = fige.add_subplot(111)    
            s = list(species).index(spec)
            axe.plot(-species_age[s], speciesbars[s], color = cladecolor[q], marker = 'o', label = spec)
            axe.legend(prop={'size':5.5})
            axe.spines['top'].set_visible(False)
            axe.spines['right'].set_visible(False)
            axe.set_xlabel('Mya')
            axe.set_ylabel('Number specificities')
            axe.set_xticks(-species_age[s])
            axe.set_xticklabels(species_age[s])
            fige.savefig(outname+'_species_clusterevolution'+spec+'.jpg', bbox_inches = 'tight', dpi = 200)
            
    
    
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
            print "Not found"
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
            print cltype
            if len(cltype) > 1:
                cltype = cltype[:-1]
                for clt in cltype:
                    for k, kingc in enumerate(kingcluster>0):
                        if np.array_equal(kingc, clusteroverlaps[clt]):
                            print k, clustercenter[clusterid == k][0], cla, ''.join(kingc.astype(int).astype(str))
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
        
        
        print 'clusterID', '#RBR', 'REP_RNCMPT', 'Ancestor', ','.join( unking), 'top10Names'
        for c, cla in enumerate(np.unique(clustage)):
            cltype = np.where(clustage == cla)[0]
            for clt in cltype:
                for k, kingc in enumerate(kingcluster):
                    if np.array_equal(kingc>0, clusteroverlaps[clt]):
                        bestrep = np.argmin(closestdist[clusterid == k])
                        #closerep = clusterrealname[np.where(clusterid == k)[0][np.argsort(closestdist[clusterid == k])]][:10]
                        print k, int(np.sum(clusterid == k)), closestmeasured[clusterid == k][bestrep], cla, ','.join(np.around(kingc).astype(int).astype(str)), ','.join(clustername[k])
        
       
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
        print np.amax(clusterparalogs)
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
        print outname
        
        fig.savefig(outname+'.jpg', bbox_inches = 'tight', dpi = 200)
        figs.savefig(outname+'_speciesclusters.jpg', bbox_inches = 'tight', dpi = 200)
        figs2.savefig(outname+'_speciesclusterslost.jpg', bbox_inches = 'tight', dpi = 200)
        figs3.savefig(outname+'_speciesevolutionscenario.jpg', bbox_inches = 'tight', dpi = 200)
        figs4.savefig(outname+'_specieshomologsequence.jpg', bbox_inches = 'tight', dpi = 200)
    plt.show()




