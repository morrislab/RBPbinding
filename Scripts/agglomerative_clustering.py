
# run: cluster_centre, clusters, cluster_assignment, cluster_scores = threshclust(identity_cutoff, sequence_names, identity_cutoff, clustering_type)

# identity_cutoff: depending on the identity matrix 0.9, or 90
# clustering_type: can be
    # 'strong': every protein in the cluster is within the identity_cutoff to each other
import numpy as np
import sys, os
from sklearn.cluster import AgglomerativeClustering
##### Clustering #####

def clustering(dmat, method, dthresh, direction):
    if direction == 'gt':
        dthresh = 1.-dthresh
    clustalg = AgglomerativeClustering(n_clusters=None, affinity='precomputed', connectivity=None, compute_full_tree=True, linkage=method, distance_threshold=dthresh).fit(dmat)
    mdb = clustalg.labels_
    clusters, mclustlen = np.unique(mdb, return_counts = True)
    print( mclustlen)
    mclustscore = []
    mclustercomp = []
    mclustcentre = []
    for m, cl in enumerate(clusters):
        mclustscore.append(np.mean(dmat[np.triu_indices(len(dmat), 1)]))
        clustind = np.where(mdb == cl)[0]
        mclustercomp.append(clustind)
        mclustcentre.append(clustind[np.argmin(np.sum(dmat[clustind][:,clustind], axis = 1))])
    return mclustcentre, mclustercomp, mdb, mclustscore, mclustlen





#### Read in function #####

def readinclustermat(clustermat):
    if os.path.splitext(clustermat)[1] == '.npz':
        files = np.load(clustermat)
        fnames = list(files.files)
        clusteringmat = files[fnames[0]]
        clustnames = files['names']
    else:
        clustnames = open(clustermat, 'r').readline().strip().split()[1:]
        if '||' in clustnames[0]:
            for c, cname in enumerate(clustnames):
                cname = cname.split('||')
                clustnames[c] = cname[-1]
        clusteringmat = np.genfromtxt(clustermat, dtype = str)
        if clusteringmat[0,1] in  clustnames:
            clusteringmat = clusteringmat[1:]
        if clusteringmat[1,0] in  clustnames:
            clusteringmat = clusteringmat[:, 1:]
        clusteringmat = clusteringmat.astype(float)
    return np.array(clustnames), clusteringmat

def keepfunc(names, matrix, rightnames):
    keep = []
    for rightname in rightnames:
        keep.append(list(names).index(rightname))
    names = np.array(names)[keep]
    matrix = matrix[keep]
    matrix = matrix[:,keep]
    return names, matrix



if __name__ == '__main__':
    simatrix = sys.argv[1] # sim, dist: sim will be translated to distance matrix
    matdir = sys.argv[2]
    insimcut = sys.argv[3]
    clusttype = sys.argv[4]  #'ward', 'complete', 'average', 'single'
    direction = sys.argv[5] # lt, gt
    
    simcut = float(insimcut)
    print( 'Read in ...')
    simnames, simat = readinclustermat(simatrix)
    simat = np.array(simat)
    if matdir == 'sim':
        simat[np.diag_indices(len(simat))] = np.amax(simat)
        simcut = simcut/np.amax(simat)
        simat = simat/np.amax(simat)
        simat = 1. - simat
        print('Reverse simcut:', simcut)

    print( 'done.')
    
    if '--proteinlist' in sys.argv:
        proteinlist = sys.argv[sys.argv.index('--proteinlist')+1]
        prots = np.sort(np.genfromtxt(proteinlist, dtype = str))
        if len(np.shape(prots)) == 2:
            prots = prots[:,0]
        simnames, simat = keepfunc(simnames, simat, prots)
    
    if '--clusters' in sys.argv:
        clfile = sys.argv[sys.argv.index('--clusters')+1]
        outname = os.path.splitext(simatrix)[0]+'_aggclusters-'+clusttype+direction+insimcut+'_on-'+os.path.splitext(os.path.split(clfile)[1])[0]+'.cl'
        clobj = np.load(clfile)
        cnnames = clobj['names']
        inclusters = clobj['clusters']
        uncluster = np.unique(inclusters)
        mclustcentre= []
        mclustercomp = []
        mdb = -np.ones(len(simnames), dtype = int)
        mclustlen = []
        c = 0
        for i, inc in enumerate(uncluster):
            icnames = cnnames[inclusters == inc]
            lsimloc = np.where(np.isin(simnames, icnames))[0]
            print(i, lsimloc)
            if len(lsimloc) > 1:
                lsimmat = simat[lsimloc][:, lsimloc]
                msclustcentre, msclustercomp, msdb, msclustscore, msclustlen = clustering(lsimmat, clusttype, simcut, direction)
                for m in range(len(msclustercomp)):
                    print (lsimloc[msclustercomp[m]])
                    mclustcentre.append(lsimloc[msclustcentre[m]])
                    mclustercomp.append(lsimloc[msclustercomp[m]])
                    mdb[lsimloc[msclustercomp[m]]] = c
                    c += 1
                    mclustlen.append(len(msclustercomp[m]))
            else:
                mclustcentre.append(lsimloc[0])
                mclustercomp.append(lsimloc)
                mdb[lsimloc[0]] = c
                mclustlen.append(1)
                c += 1 
    else:
        outname = os.path.splitext(simatrix)[0]+'_aggclusters-'+clusttype+direction+insimcut+'.cl'
        mclustcentre, mclustercomp, mdb, mclustscore, mclustlen = clustering(simat, clusttype, simcut, direction)
    
    
    print(len(mclustcentre), 'clusters found for', len(mdb), 'objects with', clusttype, insimcut, mdb)
    clen, clenum = np.unique(mclustlen, return_counts = True)
    for c in range(len(clen)):
        print(clen[c], clenum[c])
    
    if '--print_clusters' in sys.argv:
        namefile = sys.argv[sys.argv.index('--print_clusters')+1]
        realnames = np.genfromtxt(namefile, dtype = str)
        for c, cl in enumerate(mclustercomp):
            print( c, ':', len(cl))
            print( ' '.join(realnames[np.isin(realnames[:,1], simnames[cl]),2]))
        
    
    np.savez_compressed(os.path.splitext(outname)[0]+'.npz', names = simnames, clusters = mdb, clustercentre = mclustcentre, clustlist = mclustercomp)
   
