import numpy as np
import sys,os
import matplotlib.pyplot as plt
from functools import reduce
from sklearn.cluster import DBSCAN
import matplotlib.cm as cm
import matplotlib as mpl
from sklearn import preprocessing
from scipy.cluster.hierarchy import dendrogram, linkage
import sklearn.decomposition as skd
import sklearn.manifold as skm
import copy as cp
import umap
from scipy.linalg import svd
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors




def keepfunc(names, matrix, rightnames):
    keep = []
    klocation = []
    for ri, rightname in enumerate(rightnames):
        if rightname in names:
            keep.append(list(names).index(rightname))
            klocation.append(ri)
    names = np.array(names)[keep]
    matrix = matrix[keep]
    matrix = matrix[:,keep]
    return names, matrix, klocation

def normdata(clmat):
    if '--Minmaxnorm' in sys.argv:
        min_max_scaler = preprocessing.MinMaxScaler()
        clmat = min_max_scaler.fit_transform(clmat)
    else:
        clmat = preprocessing.scale(clmat)
    return clmat

##### generate coordinates for plotting with chosen embedding
def reduce2d(pfeature, pfmethod, pfdim, pfparms):
    if pfmethod == 'None':
        if np.isscalar(pfdim):
            orep2d = pfeature[:, :pfdim]
        else:
            orep2d = pfeature[:, pfdim]
    elif pfmethod == "tsne":
        if pfparms is not None:
            perpxty = float(pfparms[0])
            metic = pfparms[1]
            early_exaggeration = float(pfparms[2])
            randstat = int(pfparms[3])
        else:
            perpxty = 50.
            metic = 'euclidean'
            early_exaggeration = 12.
            randstat = 42
        print "TSNE ... "
        tsne = skm.TSNE(n_components=pfdim, perplexity=perpxty, early_exaggeration=early_exaggeration, learning_rate=200.0, n_iter=1000, n_iter_without_progress=300, min_grad_norm=1e-07, metric=metic, init='random', verbose=0, random_state=randstat, method='barnes_hut', angle=0.5)
        if len(pfeature) == 2:
            tsne.fit(pfeature[0])
            orep2d = tsne.fit_transform(pfeature[0])
            orep2db = tsne.fit_transform(pfeature[1])
            orep2d = np.append(orep2d, orep2db, axis = 0)
        else:
            orep2d = tsne.fit_transform(pfeature)
        
        
    elif pfmethod == 'svd':
        print 'SVD'
        orep2d, sig2d, fac2d = svd(pfeature, full_matrices=False)
        orep2d = orep2d[:, : pfdim]*sig2d[:pfdim]
    elif pfmethod == 'pca':
        print 'PCA ... '
        pca = skd.PCA(n_components=pfdim, copy=True, whiten=False, svd_solver='auto', tol=0.0, iterated_power='auto', random_state=None)
        orep2d = pca.fit_transform(pfeature)
    elif pfmethod == 'kernelpca':
        print 'kernelPCA ... '
        if pfparms is None:
            pcakern = 'linear'
        else:
            pcakern = pfparms # linear | poly | rbf | sigmoid | cosine | precomputed
        pca = skd.KernelPCA(n_components=pfdim, kernel=pcakern, gamma=None, degree=3, coef0=1, kernel_params=None, alpha=1.0, fit_inverse_transform=False, eigen_solver='auto', tol=0, max_iter=None, remove_zero_eig=False, random_state=None, copy_X=True, n_jobs=1)
        orep2d = pca.fit_transform(pfeature)
    elif pfmethod == 'sparsepca':
        print 'sparcePCA ... '
        if pfparms is None:
            alph = 1
        else:
            alph = float(pfparms)
        pca = skd.SparsePCA(n_components=pfdim, alpha=alph, ridge_alpha=0.01, max_iter=1000, tol=1e-08, method='lars', n_jobs=1, U_init=None, V_init=None, verbose=False, random_state=None)
        orep2d = pca.fit_transform(pfeature)
    elif pfmethod == 'truncatedpca':
        print 'truncatedPCA ... '
        pca = skd.TruncatedSVD(n_components=pfdim, algorithm='randomized', n_iter=5, random_state=None, tol=0.0)
        orep2d = pca.fit_transform(pfeature)
    elif pfmethod == 'nmf':
        if pfparms is None:
            alph = 0.
        else:
            alph = float(pfparms)        
        print 'NMF ... '
        pfeature = preprocessing.MinMaxScaler().fit_transform(pfeature)
        pca = skd.NMF(n_components=pfdim, init=None, solver='cd', beta_loss='frobenius', tol=0.0001, max_iter=200, random_state=None, alpha=alph, l1_ratio=0.0, verbose=0, shuffle=False)
        orep2d = pca.fit_transform(pfeature)
    elif pfmethod == 'isomap': 
        if pfparms is None:
            alph = 5
        else:
            alph = int(pfparms) 
        pca = skm.Isomap(n_neighbors=5, n_components=pfdim, eigen_solver='auto', tol=0, max_iter=None, path_method='auto', neighbors_algorithm='auto', n_jobs=1)
        orep2d = pca.fit_transform(pfeature)
    elif pfmethod == 'mds': # multidimensional scaling
        if pfparms is None:
            alph = True
        else:
            alph = False 
        pca = skm.MDS(n_components=pfdim, metric=alph, n_init=4, max_iter=300, verbose=0, eps=0.001, n_jobs=4, random_state=None, dissimilarity='euclidean')
        orep2d = pca.fit_transform(pfeature)
    elif pfmethod == 'umap':
        if pfparms is None:
            neigh = 14
            met = 'cosine'
            mindist = 0.7
            spread = 1.
            tseed = 42
        else:
            neigh = int(pfparms[0])
            met = pfparms[1]# euclidean, cosine, etc.
            mindist = float(pfparms[2]) 
            spread = float(pfparms[3]) 
            tseed = int(pfparms[4])
        pca = umap.UMAP(n_neighbors=neigh, n_components=pfdim, metric=met, n_epochs=None, learning_rate=1.0, init='spectral', min_dist=mindist, spread=spread, set_op_mix_ratio=1.0, local_connectivity=1.0, repulsion_strength=1.0, negative_sample_rate=5, transform_queue_size=4.0, a=None, b=None, random_state=2, metric_kwds=None, angular_rp_forest=False, target_n_neighbors=-1, target_metric='categorical', target_metric_kwds=None, target_weight=0.5, transform_seed=tseed, verbose=False)
        if len(pfeature) == 2:
            pca.fit(pfeature[0])
            orep2d = pca.transform(pfeature[0])
            orep2db = pca.transform(pfeature[1])
            orep2d = np.append(orep2d, orep2db, axis = 0)
        else:
            orep2d = pca.fit_transform(pfeature)
        

    return orep2d





nucleotideabbreviations = ['A','C','G','U','M','R','W','S','Y','K','V','H','D','B','N']
nucabbmat = np.array([[1,0,0,0],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,1],
            [1,1,0,0],
            [1,0,1,0],
            [1,0,0,1],
            [0,1,1,0],
            [0,1,0,1],
            [0,0,1,1],
            [1,1,1,0],
            [1,1,0,1],
            [1,0,1,1],
            [0,1,1,1],
            [0,0,0,0]])

scormotmat = np.copy(nucabbmat)
scormotmat[-1,:]=1.
scormotmat = (scormotmat.T.astype(float)/np.sum(scormotmat.astype(float), axis = 1)).T
# make IUPAC abbreviation of pwm
def abbreviatepwm(pwm):
    motabb = ''
    for pwnum, pwpos in enumerate(pwm):
        sortpw = np.argsort(-pwpos)
        abbvec = np.zeros(4)
        if pwpos[sortpw[0]] >= 0.47:
            if pwpos[sortpw[1]] < 0.8*0.47:
                abbvec[sortpw[0]] = 1
            else:
                abbvec[sortpw[:2]] = 1
        elif pwpos[sortpw[0]] >= 0.3:
            if pwpos[sortpw[1]] < 0.8*.3:
                abbvec[:] = 0
            elif pwpos[sortpw[1]] >= 0.8*0.3:
                if pwpos[sortpw[2]] < 0.8*0.3:
                    abbvec[sortpw[:2]] = 1
                elif pwpos[sortpw[2]] >= 0.8*0.3:
                    abbvec[sortpw[:3]] = 1
        abbrev = (nucabbmat == abbvec).all(axis = 1).nonzero()[0]
        #print abbrev, abbvec
        motabb += nucleotideabbreviations[abbrev[0]]
    return motabb

# remove N's at end of abbreviation
def cleanedgesmotif(nmotif):
    #print nmotif
    cleanl = 0
    cleanr = 0
    stend = [0,len(nmotif)-1]
    for nn in range(len(nmotif)):
        if cleanl == 0:
            if nmotif[nn] != 'N':
                cleanl += 1
        if cleanl == 1 :
            stend[0] = nn
            cleanl += 1
        if cleanr == 0:
            if nmotif[-nn-1] != 'N':
                cleanr += 1
        if cleanr == 1 :
            stend[1] = len(nmotif)-nn
            cleanr += 1
    #print stend
    nmotif = nmotif[stend[0]:stend[1]]
    #print nmotif
    return nmotif











# Generate common pwm from set of pwms
def clustermot(pwms, idnames, mclustercomp):
    clpwms = []
    clmots = []
    for c, clust in enumerate(mclustercomp):
        cpwm = combinepwms(pwms, clust)
        cmot = abbreviatepwm(cpwm)
        cmot = cleanedgesmotif(cmot)
        clpwms.append(cpwm)
        clmots.append(cmot)
    return clpwms, clmots

def kldiv(a1, a2):
    out = np.sum(a1*np.log(a1/a2)+a2*np.log(a2/a1))
    return out

def combinepwms(pset, clust):
    psetclust = []
    for clu in clust:
        psetclust.append(pset[clu])
    while len(psetclust) > 1:
        klds = []
        klpwms = []
        pairs = []
        for c, clu in enumerate(psetclust):
            for d in range(c+1, len(psetclust)):
                clw = psetclust[d]
                kld, klpwm = alignpwms(clu, clw)
                pairs.append([c,d])
                klds.append(kld)
                klpwms.append(klpwm)
        bekl = np.argmin(klds)
        npsclust = []
        for c, clu in enumerate(psetclust):
            if c not in pairs[bekl]:
                npsclust.append(clu)
        npsclust.append(klpwms[bekl])
        psetclust = npsclust
    return psetclust[0]
            
    
def meanpwm(pa, pb, ofs):
    la = len(pa)
    lb = len(pb)
    pfin = np.zeros((len(pb),4))
    pfin += pb
    pfin[min(0,-(lb+ofs)):min(lb,la-ofs)] += pa[max(0,ofs):min(la,ofs+lb)]
    pnorm = np.ones(lb)
    pnorm[min(0,-(lb+ofs)):min(lb,la-ofs)] +=1.
    return (pfin.T/pnorm).T

def alignpwms(x, y):
    lenx = len(x)
    leny = len(y)
    if lenx > leny:
        a = x
        b = y
    else:
        a = y
        b = x
    #print a
    #print b
    klset = []
    offset = []
    for i in range(len(a)-len(b)+1):
        pcheck1 = np.copy(b)
        pcheck2 = np.ones((len(b), 4))*0.25
        pcheck2 = a[i:i+len(b)]
        klset.append(kldiv(pcheck1, pcheck2))
        offset.append(i)
    for i in range(2,len(b)):
        pcheck1 = np.copy(b)
        pcheck2 = np.ones((len(b), 4))*0.25
        pcheck2[-i:] = a[:i]
        klset.append(kldiv(pcheck1, pcheck2))
        offset.append(i-len(b))
        pcheck3 = np.ones((len(b), 4))*0.25
        pcheck3[:i] = a[-i:]
        klset.append(kldiv(pcheck1, pcheck3))
        offset.append(len(a)-i)
    bl = np.argmin(klset)
    bestoff = offset[bl]
    bestkl = klset[bl]
    mpwm = meanpwm(a, b, bestoff)
    return bestkl, mpwm





#Coloring according to motif
def motcolor(motiflist, motctype):
    #import matplotlib.colors as allcolors
    #colorlist = allcolors.cnames
    colormat = []
    dissolve = np.ones(11)
    dissolvemix = np.ones(4)
    if motctype == 'piechart' or 'manual':
        for m, cmotif in enumerate(motiflist):
            cl = 0
            ntdist = np.zeros(4)
            for cmot in cmotif:
                if cmot != 'N':
                    cl += 1.
                if cmot == 'A':
                    # green
                    ntdist[0]+=1.
                if cmot == 'C':
                    # blue
                    ntdist[1]+=1.
                if cmot == 'G':
                    # yellow
                    ntdist[2]+=1.
                if cmot == 'U':
                    # red
                    ntdist[3]+=1.
                if cmot == 'M':
                    ntdist[0] += 0.5
                    ntdist[1] += 0.5
                if cmot == 'K':
                    ntdist[2] += 0.5
                    ntdist[3] += 0.5
                if cmot == 'W':
                    ntdist[0] += 0.5
                    ntdist[3] += 0.5
                if cmot == 'R':
                    ntdist[0] += 0.5
                    ntdist[2] += 0.5
                if cmot == 'Y':
                    ntdist[1] += 0.5
                    ntdist[3] += 0.5
                if cmot == 'B':
                    ntdist[1] += 0.33
                    ntdist[2] += 0.33
                    ntdist[3] += 0.33
                if cmot == 'S':
                    ntdist[1] += 0.5
                    ntdist[2] += 0.5
                if cmot == 'D':
                    ntdist[0] += 0.33
                    ntdist[2] += 0.33
                    ntdist[3] += 0.33
                if cmot == 'H':
                    ntdist[1] += 0.33
                    ntdist[0] += 0.33
                    ntdist[3] += 0.33
                if cmot == 'V':
                    ntdist[1] += 0.33
                    ntdist[2] += 0.33
                    ntdist[0] += 0.33
            
            if cl != 0 and motctype == 'piechart':
                colormat.append(np.around(ntdist/cl, 2))
            elif cl != 0 and motctype == 'manual':
                colormat.append(ntdist)
            else:
                colormat.append(np.ones(4)/4.)
                
        if motctype == 'manual':
            #print np.argsort(np.argsort(colormat))
            tab40 = np.concatenate([cm.tab20b([0,2]),cm.tab20c(np.arange(0,2, dtype = int)), cm.tab20([18]), cm.tab20c(np.arange(8,10, dtype = int)),cm.tab20b(np.arange(4,6, dtype = int)),cm.tab20([16]), cm.tab20b(np.arange(8,11, dtype = int)),cm.tab20([6,10]),cm.tab20c([4,5]),cm.Dark2([3]),cm.tab20b(np.arange(16,19, dtype = int)), cm.gist_rainbow(np.linspace(0,1,10))], axis = 0)
            colormat = tab40[np.argsort(np.argsort(np.array(colormat)[:,3]+np.array(colormat)[:,1]))]
            
    colormat = np.array(colormat)
    return colormat




#figsettings= [fwidth, fheight, falpha, fscatsize, axisvisual, fconsize, fhead, anfontsize]
def plotevo(rep2d, annotatemot, annotate, annotateloc, scattercolor, colorscheme, clustcentre, clustercomp, clustervisual,shape, connectmap, protext, plotname, figsettings):
    ### open figure
    
    edgeclist = ['k', 'none', 'none', 'none']
    
    if plotname is not None:
        figurename = os.path.splitext(plotname)[0]
        figfmt = os.path.splitext(plotname)[-1]
    else:
        figurename = 'Motif scatter'
    
    fig = plt.figure(figurename, figsize=(figsettings[0],figsettings[1]), dpi = 100)
    dimensions = len(rep2d[0])
    axlim = [np.amin(rep2d, axis = 0), np.amax(rep2d, axis = 0)]
    diaglen = (axlim[1]- axlim[0])**2
    centre = np.mean(rep2d, axis = 0)
    
    
    if clustervisual != 'None':
        clusterlen = []
        for comp in enumerate(clustercomp):
            clusterlen.append(len(comp))
        clusterlen = np.array(clusterlen)
    
    if annotateloc == 'outside':
        figanno = plt.figure(figurename+'annotation', figsize=(figsettings[0]*.4,figsettings[1]), dpi = 100)
        axanno = figanno.add_subplot(111)
        axanno.spines['left'].set_visible(False)
        axanno.spines['top'].set_visible(False)
        axanno.spines['bottom'].set_visible(False)
        axanno.spines['right'].set_visible(False)
        axanno.tick_params( axis='x',which='both', bottom=False, top=False, labelbottom=False, labeltop = False)
        axanno.tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
        yanno = np.sum(annotate)
        axanno.set_xlim([-1.,2.])
        ycount = 0
        for ii, manno in enumerate(annotatemot):
            if annotate[ii]:
                if colorscheme != 'piechart':
                    axanno.scatter([0], yanno-ycount, c = scattercolor[ycount], alpha = figsettings[2], s = figsettings[3]*150)
                else:
                    colorscatter([0], yanno-ycount, scattercolor[ycount], figsettings[3]*0.8, axanno, 0.3, figsettings[3])
                axanno.annotate(manno + ' '+str(len(clustercomp[ii])), (0.5, yanno-ycount), horizontalalignment='left', verticalalignment='center', fontsize = figsettings[-1], fontweight = None) #'bold' )
                ycount +=1
        
    for d2 in range(dimensions):
        for d1 in range(d2+1, dimensions):
            ax = fig.add_subplot(dimensions-1, dimensions-1, d2*(dimensions-1) + d1)
            #ax = fig.add_subplot2grid(gridim, (d1, d2))
            ax.scatter(centre[d1], centre[d2], color = 'k', marker = 'x')
            if figsettings[4] == False:
                ax.tick_params( axis='x',which='both', bottom=False, top=False, labelbottom=False, labeltop = False)
                ax.tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
                ax.spines['left'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['right'].set_visible(False)
            
            if len(connectmat) > 0:
                maxlen = np.amax(np.sqrt(rep2d[:,d1]**2 + rep2d[:,d2]**2))/.5
                for r in range(len(connectmat)):
                    for s in range(r+1, len(connectmat)):
                        distrs = np.sqrt((rep2d[r, d1] - rep2d[s, d1])**2 + (rep2d[r, d2]- rep2d[s, d2])**2)
                        if (distrs < maxlen) and connectmat[r,s] > 0:
                            ax.plot([rep2d[r, d1], rep2d[s, d1]], [rep2d[r, d2], rep2d[s, d2]], c = 'grey', alpha = connectmat[r,s]/(1.5*np.amax(connectmat)), lw = 2.*connectmat[r,s]/np.amax(connectmat))
                
            
            if clustervisual == 'None':
                annoloc = rep2d[:, [d1,d2]]
                if colorscheme != 'piechart':
                    ax.scatter(rep2d[:,d1], rep2d[:,d2], c= scattercolor, alpha = figsettings[2], s = figsettings[3])
                else:
                    for i in range(len(rep2d)):
                        colorscatter(rep2d[i,d1], rep2d[i,d2], scattercolor[i], figsettings[3], ax, figsettings[2], figsettings[3])
            
            
            elif clustervisual == 'centresize':
                annoloc = rep2d[clustcentre][:, [d1,d2]]
                if colorscheme != 'piechart':
                    ax.scatter(rep2d[clustcentre,d1], rep2d[clustcentre,d2], c = scattercolor, alpha = figsettings[2], s = figsettings[3]*clusterlen)
                else:
                    for i in range(len(clustcentre)):
                        colorscatter(rep2d[clustcentre[i],d1], rep2d[clustcentre[i],d2], scattercolor[i], clusterlen[i], ax, figsettings[2], figsettings[3])
                
            elif clustervisual == 'centrecolor':
                annoloc = rep2d[clustcentre][:, [d1,d2]]
                restcluster = np.concatenate(clustercomp)
                notclustcentre = np.delete(np.arange(len(rep2d)), restcluster)
                #### add coloring of background points 
                for sh in np.unique(shape):
                    shapemask = shape[notclustcentre] == sh
                    ax.scatter(rep2d[notclustcentre,d1][shapemask], rep2d[notclustcentre,d2][shapemask], c = 'grey', edgecolor = 'none', alpha = 0.3, s = figsettings[3]*90, marker = markerlist[sh])
                for c, ccomp in enumerate(clustercomp):
                    ccomp = ccomp[ccomp != clustcentre[c]]
                    for i in ccomp:
                        if colorscheme != 'piechart':
                            ax.scatter(rep2d[i,d1], rep2d[i,d2], c = scattercolor[c], alpha = figsettings[2], s = figsettings[3]*120, marker = markerlist[shape[i]], edgecolor = edgeclist[shape[i]])
                        else:
                            colorscatter(rep2d[i,d1], rep2d[i,d2], scattercolor[c], figsettings[3]*0.8, ax, 0.3, figsettings[3])
                        #ax.scatter(rep2d[ccomp,d1], rep2d[ccomp,d2], c = scattercolor[c] , alpha = 0.2, s = figsettings[3]*100)
                ### Need to reestablish for some stuff
                if colorscheme != 'piechart':
                    for sh in np.unique(shape):
                        shapemask = shape[clustcentre] == sh
                        clustcentre = np.array(clustcentre)
                        ax.scatter(rep2d[clustcentre[shapemask],d1], rep2d[clustcentre[shapemask],d2], c = scattercolor[shapemask], alpha = figsettings[2], s = figsettings[3]*140., marker = markerlist[sh], edgecolor = edgeclist[sh])
                else:
                    for i in range(len(clustcentre)):
                        colorscatter(rep2d[clustcentre[i],d1], rep2d[clustcentre[i],d2], scattercolor[i], figsettings[3]*3., ax, figsettings[2], figsettings[3]*multfact)
                
            elif clustervisual == 'meansize':
                meanrep2d = []
                for comp in clustercomp:
                    meanrep2d.append(np.mean(rep2d[comp], axis = 0))
                meanrep2d = np.array(meanrep2d)
                annoloc = meanrep2d[:, [d1,d2]]
                if colorscheme != 'piechart':
                    ax.scatter(meanrep2d[:,d1], meanrep2d[:,d2], c = scattercolor, alpha = figsettings[2], s = figsettings[3]*clusterlen)
                else:
                    for i in range(len(clustcentre)):
                        colorscatter(meanrep2d[i,0], meanrep2d[i,1], scattercolor[i], clusterlen[i], ax, figsettings[2], figsettings[3])
                
            elif clustervisual == 'meancolor':
                meanrep2d = []
                for comp in clustercomp:
                    meanrep2d.append(np.mean(rep2d[comp], axis = 0))
                meanrep2d = np.array(meanrep2d)
                annoloc = meanrep2d[:, [d1,d2]]
                notclustcentre = np.delete(np.arange(len(rep2d)), clustcentre)
                ax.scatter(rep2d[notclustcentre,d1], rep2d[notclustcentre,d2], c = 'grey', alpha = 0.2, s = figsettings[3])
                if colorscheme != 'piechart':
                    ax.scatter(meanrep2d[:,d1], meanrep2d[:,d2], c = scattercolor, alpha = figsettings[2], s = figsettings[3])
                else:
                    for i in range(len(clustcentre)):
                        colorscatter(meanrep2d[i,0], meanrep2d[i,1], scattercolor[i], 2., ax, figsettings[2], figsettings[3])
                    
            if len(protext) > 0:
                for pte in protext:
                    ax.text(pte[1][d1],pte[1][d2],pte[0], ha = 'left', va = 'bottom')
            
            
            
            if annotatemot is not None and annotateloc == 'inside':
                # clean annotation if they are too close together, could be added!!!
                for ii, manno in enumerate(annotatemot):
                    if annotate[ii]:
                        ax.annotate(manno, (annoloc[ii][0], annoloc[ii][1]), horizontalalignment='right', verticalalignment='top', fontsize = figsettings[-1], fontweight = None) #'bold' )
    
    if '--savefig' in sys.argv:
        
        fig.savefig(figurename+figfmt, format=figfmt.strip('.'), bbox_inches='tight', dpi = 300)
        if annotateloc == 'outside':
            figanno.savefig(figurename+'_annotation'+figfmt, format=figfmt.strip('.'), dpi = 200)
    else:
        plt.show()

# color by eigenvector strength
def plotfactor(rep2d, colorfactor, numfactors, combinedfactor, nameprots, normalize, plotname, figsettings):
    
    calpha = 1.
    if combinedfactor:
        fig = plt.figure(plotname, figsize=(figsettings[0]*.5,figsettings[1]*0.5), dpi = 100)
        ax = []
    dimensions = len(rep2d[0])
    centre = np.mean(rep2d, axis = 0)
    
    nfactors = len(numfactors)
    
    if plotname is not None:
        figurename = os.path.splitext(plotname)[0]+'-factorcolor'
        figfmt = os.path.splitext(plotname)[-1]
    else:
        figurename = 'Latent_scatter_factorcolor'
    
    for co, colorfac in enumerate(colorfactor):
        if normalize:
            colorfac[colorfac<0] = colorfac[colorfac<0]/np.amax(np.absolute(colorfac[colorfac<0]))
            colorfac[colorfac>0] = colorfac[colorfac>0]/np.amax(np.absolute(colorfac[colorfac>0]))
        # generate colormap
        if combinedfactor:
            lowcolor = np.array([[0.6,0.1,0.6,1.], [0.8,0.45,0.1,1.], [0.2,0.6,0.6,1.],[0.2,0.2,0.6,1.],[0.6,0.6,0.2,1.],[0.6,0.2,0.2,1.],[0.3,0.8,0.7,1.],[0.1,0.8,0.4,1.]])[co]
            lowcolor = np.array(lowcolor)
            highcolor = np.ones(4) - lowcolor
            lowcolor = np.array(lowcolor)
            highcolor = np.array(highcolor)
            darkness = 0.7
            lowcolor = lowcolor*darkness/np.amax(lowcolor)
            highcolor = highcolor*darkness/np.amax(highcolor)
            highcolor[3] = 1.
            lowcolor[3] = 1.
            

        else:
            lowcolor = np.array([0.45,0.1,0.6,1.])
            highcolor = np.array([0.2,0.6,0.2,1.])
        
        if (co == 0 and combinedfactor) or not combinedfactor:
            middlecolor = np.array([0.9,0.9,0.9,1.])
        else:
            middlecolor = np.array([0.95,0.95,0.95,0.])
        mixcolorlow = (3.*middlecolor + lowcolor)/4.
        mixcolorhigh = (3.*middlecolor + highcolor)/4.
        mixcolorlowb = (1.*middlecolor + 3.*lowcolor)/4.
        mixcolorhighb = (1.*middlecolor + 3.*highcolor)/4.
        tricol = [lowcolor, mixcolorlowb, mixcolorlow, middlecolor, mixcolorhigh, mixcolorhighb,highcolor]
        cvals = [-np.amax(np.absolute(colorfac)),-2*np.amax(np.absolute(colorfac))/3.,-1.*np.amax(np.absolute(colorfac))/3.,0., 1*np.amax(np.absolute(colorfac))/3.,  2*np.amax(np.absolute(colorfac))/3., np.amax(np.absolute(colorfac))]
        norm=plt.Normalize(min(cvals),max(cvals))
        tuples = list(zip(map(norm,cvals), tricol))
        cmapfactor = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples) 
        
        
        
        if not combinedfactor:
            if plotname is not None:
                figurename = os.path.splitext(plotname)[0]+'-factorcolor'+str(numfactors[co])
                figfmt = os.path.splitext(plotname)[-1]
            else:
                figurename = 'Latent_scatter_factorcolor'+str(numfactors[co])
            fig = plt.figure(figurename, figsize=(figsettings[0]*.5,figsettings[1]*0.5), dpi = 100)
            ax = []
        else:
            figurename += '-'+str(numfactors[co])
        vmin = -np.amax(np.absolute(colorfac))
        vmax = np.amax(np.absolute(colorfac))
            
        if combinedfactor:
            colorbar = fig.add_subplot(nfactors, 1, co+1)
            colorbar.set_position([0.33,0.-0.055*co,0.33,0.05])
            colorbar.tick_params(left = False, labelleft = False, labelbottom = True, bottom = True, right = False, labelright = True)
            colorbar.set_yticks([.0])
            colorbar.set_yticklabels(['Z'+str(numfactors[co])])
            if co == nfactors-1:
                colorbar.set_xticks([-.5,99.5])
                colorbar.set_xticklabels(['Min', 'Max'])
                
            else:
                colorbar.tick_params(left = False, labelleft = False, labelbottom = False, bottom = False, right = False)
        else:
            colorbar = fig.add_subplot(nfactors, 1, co+1)
            colorbar.set_position([0.33,0.,0.33,0.05])
            colorbar.set_xticks([-0.5, 99.5])
            colorbar.set_xticklabels(['Min', 'Max'])
            colorbar.tick_params(left = False, labelleft = False, right = False, labelright = False)
        colorbar.imshow([np.linspace(vmin,vmax, 100)], origin = 'lower', cmap = cmapfactor, aspect = 'auto', norm = norm, vmin = vmin, vmax = vmax)
        
        
        sx = 0 
        for d2 in range(dimensions):
            for d1 in range(d2+1, dimensions):
                if (combinedfactor and co == 0) or not combinedfactor:
                    ax.append(fig.add_subplot(dimensions-1, dimensions-1, d2*(dimensions-1) + d1))
                    if figsettings[4] == False:
                        ax[sx].tick_params( axis='x',which='both', bottom=False, top=False, labelbottom=False, labeltop = False)
                        ax[sx].tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
                        ax[sx].spines['left'].set_visible(False)
                        ax[sx].spines['top'].set_visible(False)
                        ax[sx].spines['bottom'].set_visible(False)
                        ax[sx].spines['right'].set_visible(False)
                
                colorsort = np.argsort(np.absolute(colorfac))
                p1 = ax[sx].scatter(rep2d[:,d1][colorsort], rep2d[:,d2][colorsort], c= colorfac[colorsort] , cmap = cmapfactor, norm = norm, s = figsettings[3]*80)
                ax[sx].scatter(centre[d1], centre[d2], color = 'k', marker = 'x')
                if nameprots is not None:
                    for pnrt in nameprots[co]:
                        ax[sx].text(pnrt[0][d1], pnrt[0][d2], pnrt[1])
                sx += 1
        if not combinedfactor or (combinedfactor and co == nfactors-1):
            if '--savefig' in sys.argv:
                print figurename+figfmt
                fig.savefig(figurename+figfmt, format=figfmt.strip('.'), bbox_inches = 'tight', dpi = 300)
            else:
                plt.show()
        plt.close()
    


# color in 4 blocks according to fraction of A,C,G,U
def colorscatter(x0, x1, coloring, size, subplt, scatalph, scatsize):
# first define the ratios, ratios of colors are ap, cp, gp, up
    ap, cp, gp, up = coloring
    aT = False
    cT = False
    gT = False
    uT = False
    if ap != 0.:
        aT = True
    if cp != 0.:
        cT = True
    if gp != 0.:
        gT = True
    if up != 0.:
        uT = True
    ap = ap
    cp = ap + cp
    gp = cp + gp
    up = 1.
    # define some sizes of the scatter marker
    sizes = np.array([80, 100, 140]) * size * scatsize
    # calculate the points of the first pie marker
    # these are just the origin (0,0) +
    # some points on a circle cos,sin
    x = [0] + np.cos(np.linspace(0, 2 * np.pi * ap, 100)).tolist()
    y = [0] + np.sin(np.linspace(0, 2 * np.pi * ap, 100)).tolist()
    xy1 = np.column_stack([x, y])

    x = [0] + np.cos(np.linspace(2 * np.pi * ap, 2 * np.pi * cp, 100)).tolist()
    y = [0] + np.sin(np.linspace(2 * np.pi * ap, 2 * np.pi * cp, 100)).tolist()
    xy2 = np.column_stack([x, y])

    x = [0] + np.cos(np.linspace(2 * np.pi * cp, 2 * np.pi * gp, 100)).tolist()
    y = [0] + np.sin(np.linspace(2 * np.pi * cp, 2 * np.pi * gp ,100)).tolist()
    xy3 = np.column_stack([x, y])
    
    x = [0] + np.cos(np.linspace(2 * np.pi * gp, 2 * np.pi, 100)).tolist()
    y = [0] + np.sin(np.linspace(2 * np.pi * gp, 2 * np.pi ,100)).tolist()
    xy4 = np.column_stack([x, y])
    
    #print ap, cp, gp, up
    #print aT, cT, gT, uT
    #print x0,x1
    #print xy1, xy2, xy3, xy4
    if aT:
        subplt.scatter([x0], [x1], marker=xy1,
            s=sizes, facecolor='green', edgecolor = 'none', alpha = scatalph)
    if cT:
        subplt.scatter([x0], [x1], marker=xy2,
            s=sizes, facecolor='blue', edgecolor = 'none',alpha = scatalph)
    if gT:
        subplt.scatter([x0], [x1], marker=xy3,
            s=sizes, facecolor='yellow', edgecolor = 'none',alpha = scatalph)
    if uT:
        subplt.scatter([x0], [x1], marker=xy4,
            s=sizes, facecolor='red', edgecolor = 'none',alpha = scatalph)

    return







def readin(plotfeatures, featureaxis, prots, reduceb):
    if os.path.splitext(plotfeatures)[-1] == '.npz':
        obj = np.load(plotfeatures, allow_pickle = True)
        ofiles = obj.files
        if 'features' in ofiles:
            plotfeatures = obj['features']
            pnames = obj['expnames']
        elif 'profiles' in ofiles:
            plotfeatures = obj['profiles']
            pnames = obj['names']
            
    else:
        pnames = open(plotfeatures, 'r').readline().strip().split()[1:]
        if '||' in pnames[0]:
            for c, cname in enumerate(pnames):
                pnames[c] = cname.split('||')[-1]
        plotfeatures = np.genfromtxt(plotfeatures)
    
    if featureaxis == '0':
        plotfeatures = plotfeatures.T
    
    if prots is not None:
        if featureaxis == 'both':
            pnames, plotfeatures, klocation = keepfunc(pnames, plotfeatures, prots)
            np.fill_diagonal(plotfeatures, np.amax(plotfeatures))
        else:
            keep = []
            klocation = []
            for cna, cname in enumerate(prots):
                if cname in pnames:
                    keep.append(pnames.index(cname))
                    klocation.append(cna)
            pnames = np.array(pnames)[keep]
            plotfeatures = plotfeatures[keep]
    else:
        klocation = np.arange(len(pnames), dtype = int)
    if reduceb or prots is None:
        prots = pnames[:]
    return pnames, plotfeatures, prots, klocation
    

def readinpwm(pwfile):
    obj = open(pwfile,'r').readlines()
    pws = []
    pwname = []
    shortpws = []
    for l, line in enumerate(obj):
        line = line.strip().split()
        if len(line) != 0:
            if line[0] == 'Motif':
                if len(pwname) > 0:
                    pws.append(np.array(pw, dtype = float))
                pwname.append(line[1])
                pw = []
                shortpws.append(obj[l+1].split()[1])
            if line[0].isdigit():
                pw.append(line[1:])
    pws.append(np.array(pw, dtype = float))
    return pwname, pws, shortpws


def checkbool(string):
    if string == 'TRUE' or  string == 'true' or string == 'True':
        return True
    else:
        return False


















#### program starts
if __name__ == '__main__':

    
    prots = None
    if '--proteinlist' in sys.argv:
        proteinlist = sys.argv[sys.argv.index('--proteinlist')+1]
        prots = np.sort(np.genfromtxt(proteinlist, dtype = str)[:,1])
        
    # Determine the 2D coordinates of the proteins
    if '--plotfeatures' in sys.argv:
        print 'Read in plot features...'
        plotfeatures = sys.argv[sys.argv.index('--plotfeatures')+1] 
        featureaxis = sys.argv[sys.argv.index('--plotfeatures')+2] # 0, 1, both (both is matrix format)
        method2d = sys.argv[sys.argv.index('--plotfeatures')+3] # tsne, pca, umap, None
        nreduce = sys.argv[sys.argv.index('--plotfeatures')+4] # number of dimensions against each other ( pca plot, f.e. 1-2, 1-3, 2-3)
        
        if ',' in nreduce:
            nreduce = np.array(nreduce.split(','), dtype = int)
        else:
            nreduce = int(nreduce)
        
        pnames, plotfeatures, prots, protorder = readin(plotfeatures, featureaxis, prots, True)
        if '--transformfeatures' in sys.argv:
            transform = sys.argv[sys.argv.index('--transformfeatures')+1]
            tnfeat = int(sys.argv[sys.argv.index('--transformfeatures')+2]) # None just reduces number of orignial features, PCA can be used on k-mers first before tsne visulalization
            treduceparams = None
            plotfeatures = reduce2d(plotfeatures, transform, tnfeat, treduceparams)
            
            plotfeatures = plotfeatures[:,:topfeat]
        
        if '--normplotfeatures' in sys.argv:
            plotfeatures = (plotfeatures.T/np.sqrt(np.sum(plotfeatures**2, axis = 1))).T
        
        if '--secondfeatures' in sys.argv and '--transformfeatures' not in sys.argv and '--normplotfeatures' not in sys.argv:
            pnames2, plotfeatures2, prots2, protorder2 = readin(sys.argv[sys.argv.index('--secondfeatures')+1], featureaxis, None, True)
            pnames = np.append(pnames, pnames2)
            prots = np.append(prots, prots2)
            plotfeatures = [plotfeatures, plotfeatures2]
        
        if '--reduceparams' in sys.argv:
            reduceparams = sys.argv[sys.argv.index('--reduceparams')+1]
            if ',' in reduceparams:
                reduceparams = reduceparams.split(',')
        else:
            reduceparams = None
        # reduce the dimension of the plot features to nreduce (min 2d)
        rep2d = reduce2d(plotfeatures, method2d, nreduce, reduceparams)
        
    else:
        print 'Features to plot need to be determined'
        sys.exit()



    clustercomp = []
    clustcentre = []
    clustervisual = 'None'
    clusteroutname = ''
    if '--clusterassignment' in sys.argv:
        # cluster assignment based on domains for example
        clustervisual = sys.argv[sys.argv.index('--clusterassignment')+1]
        clusterfile = sys.argv[sys.argv.index('--clusterassignment')+2]
        if os.path.splitext(clusterfile)[-1] == '.npz':
            clnpz = np.load(clusterfile)
            clustnames = clnpz['names']
            clusterdbs = clnpz['clusters']
            keep = []
            clprotnumbers = []
            for pr, prot in enumerate(prots):
                if prot in clustnames:
                    keep.append(list(clustnames).index(prot))
                    clprotnumbers.append(pr)
                #else:
                    #print prot, 'not in clustnames'
                    #sys.exit()
            clustnames = clustnames[keep]
            clusterdbs = clusterdbs[keep]
            clusterlen = []
            clustercomp = []
            clustcentre = []
            clusters = np.unique(clusterdbs)
            for c, cna in enumerate(clusters):
                ccomp = np.where(clusterdbs == cna)[0]
                if len(ccomp) > 0:
                    clustercomp.append(ccomp)
                    clusterlen.append(len(ccomp))
                    clustcentre.append(ccomp[0])
            clusters = np.arange(len(clustercomp), dtype = int)
            clustercomp = np.array(clustercomp)
            db = clusterdbs
        elif os.path.splitext(clusterfile)[-1] == '.txt':
            cltxt = np.genfromtxt(clusterfile, dtype = str)
            clustnames = cltxt[:, 0]
            dbtype = cltxt[:, -1]
            clprotnumbers = []
            keep = []
            for pr, prot in enumerate(prots):
                if prot in clustnames:
                    keep.append(list(clustnames).index(prot))
                    clprotnumbers.append(pr)
            clustnames = clustnames[keep]
            dbtype = dbtype[keep]
            clusters, clustcentre, clusterlen = np.unique(dbtype, return_index = True, return_counts = True)
            db = -np.ones(len(clustnames), dtype = int)
            clustercomp = []
            for c, clust in enumerate(clusters):
                db[dbtype == clust] = c
                clustercomp.append(np.where(dbtype == clust)[0])
            clustercomp = np.array(clustercomp)
            
        clusteroutname = os.path.splitext(clusterfile)[0]+'-motifsetvisual'
        if -1 in clusters:
            maskc = clusters != -1
            clusters = clusters[maskc]
            clustercomp = clustercomp[maskc]
            clusterlen = clusterlen[maskc]
            
    if '--clusterassignment' in sys.argv:
        ### only keep large clusters for the rest
        if '--largecluster' in sys.argv:
            #print len(clusterlen)
            clusterlen = np.array(clusterlen)
            minclustsize = int(sys.argv[sys.argv.index('--largecluster')+1])
            print 'Keep only big clusters', minclustsize
            mask = np.array(clusterlen) >= minclustsize
            clusterlen = np.array(clusterlen)[mask]
            clustcentre = np.array(clustcentre)[mask]
            
            clustercomp = clustercomp[mask]
            for m, ma in enumerate(mask):
                if ma == False:
                    db[db == clusters[m]] = -1
                    # -1 is not a cluster, it stands for not assigned, therefore, does not have a length or centre
            clusters = clusters[mask]
            clusteroutname+= '-minsize'+str(minclustsize)
        if len(clustnames) < len(rep2d) or not np.array_equal(clprotnumbers, np.arange(len(clprotnumbers))):
            clprotnumbers = np.array(clprotnumbers)
            newdb = -np.ones(len(rep2d))
            newdb[clprotnumbers] = db
            db = newdb
            for cl, clco in enumerate(clustercomp):
                clustercomp[cl] = clprotnumbers[clco]
                clustcentre[cl] = clprotnumbers[clustcentre[cl]]
        
        if '--remove_nonclusters' in sys.argv:
            maskr = db != 0
            db = db[maskr]
            rep2d = rep2d[maskr]
            plotfeatures = plotfeatures[maskr]
            clprotnumbers = np.arange(len(db))
            for cl, clco in enumerate(clustercomp):
                clustercomp[cl] = clprotnumbers[clco]
                clustcentre[cl] = clprotnumbers[clustcentre[cl]]

    # clusters should be used as annotation
    
    
    connectmat = []
    if '--connect' in sys.argv:
        connectfile = sys.argv[sys.argv.index('--connect')+1]
        minconnect = float(sys.argv[sys.argv.index('--connect')+2])
        connectnames, connectmat, nprots, conprotnumbers = readin(connectfile, 'both', prots, False)
        connectmat[connectmat < minconnect] = 0.
        
    
    shape = None
    markerlist = np.array(['^','s', 'o', 'D', 'p', 'h', '*','X'])
    if '--motifshape' in sys.argv:
        motifshapefile = sys.argv[sys.argv.index('--motifshape')+1]
        shapecltxt = np.genfromtxt(motifshapefile, dtype = str)
        shapeclustnames = shapecltxt[:, 1]
        shapedbtype = shapecltxt[:, -1]
        keep = []
        for pr, prot in enumerate(prots):
            if prot in shapeclustnames:
                keep.append(list(shapeclustnames).index(prot))
        shapeclustnames = shapeclustnames[keep]
        shapedbtype = shapedbtype[keep]
        shapeclusters, shapeclustcentre, shapeclusterlen = np.unique(shapedbtype, return_index = True, return_counts = True)
        shapedb = -np.ones(len(shapeclustnames), dtype = int)
        shapeclustercomp = []
        for c, clust in enumerate(shapeclusters):
            shapedb[shapedbtype == clust] = c
            shapeclustercomp.append(np.where(shapedbtype == clust)[0])
        shape = shapedb
        
        #### can be added as additional vector to assign different shapes
    else:
        shape = np.ones(len(rep2d), dtype = int)*2
    
    colorscheme = 'smooth'
    if '--color_motifcomposition' in sys.argv:
        # assign color to protein motifs based on their top kmers
        cmotiffile = sys.argv[sys.argv.index('--color_motifcomposition')+1]
        cmottype = sys.argv[sys.argv.index('--color_motifcomposition')+2]
        colorscheme = sys.argv[sys.argv.index('--color_motifcomposition')+3]
        if cmottype == 'pwm':
            # readin pwms
            pwmnames, pwms, clustermotifs = readinpwm(cmotiffile)
            spwms = []
            coloringmot = []
            coloringloc = []
            for i, idname in enumerate(prots):
                if idname in pwmnames:
                    coloringloc.append(i)
                    spwms.append(pwms[pwmnames.index(idname)])
                    coloringmot.append(clustermotifs[pwmnames.index(idname)])
            coloringmot = np.array(coloringmot)
            if clustervisual != "None":
                clusterpwms, coloringmot = clustermot(spwms, prots, clustercomp)
                print 'Cluster PWM file written', clusteroutname+'-pwms.txt'
                pwmobj = open(clusteroutname+'-pwms.txt', 'w')
                for cp, clusterpwm in enumerate(clusterpwms):
                    pwmobj.write('Motif\tCluster_'+str(cp)+'\n#\t'+coloringmot[cp]+' '+str(len(clustercomp[cp]))+' '+','.join(prots[clustercomp[cp]])+'\n'+'Pos\tA\tC\tG\tU\n')
                    for i in range(len(clusterpwm)):
                        pwmobj.write(str(i+1)+'\t'+'\t'.join(np.around(clusterpwm[i],6).astype(str))+'\t\n')
                    pwmobj.write('\n\n')
                pwmobj.close()
                

        if colorscheme == 'smooth':
            scattercolor = motcolor(coloringmot, 'smooth')
        elif colorscheme == 'semismooth':
            print 'semismooth'
            scattercolor = motcolor(coloringmot, 'semismooth')
        elif colorscheme == 'manual':
            scattercolor = motcolor(coloringmot, 'manual')
        else:
            scattercolor = motcolor(coloringmot, 'piechart')
        
    elif '--color_individual' in sys.argv:
        scattercolor = np.arange(len(prots))
    elif '--color_cluster' in sys.argv:
        cmcolor = sys.argv[sys.argv.index('--color_cluster')+1]
        if clustervisual != "None":
            if cmcolor == 'rainbow':
                scattercolor = cm.rainbow(np.linspace(0., 1., len(clustcentre)))
            elif cmcolor == 'tab20b':
                scattercolor = cm.tab20b(np.arange(0, len(clustcentre)))
            elif cmcolor == 'manual':
                #tricol = np.array(['firebrick', 'saddlebrown', 'orange', 'gold', 'darkturquoise', 'steelblue', 'midnightblue', 'purple'])
                tricol = np.array(['silver', 'limegreen', 'green', 'grey'])
                cvals = np.unique(db)
                norm=plt.Normalize(min(cvals),max(cvals))
                tuples = list(zip(map(norm,cvals), tricol))
                cmapmanual = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)
                print tuples
                scattercolor = cmapmanual(np.linspace(0, 1, len(clustcentre)))
        else:
            if cmcolor == 'rainbow':
                scattercolor = cm.rainbow((db - np.amin(db))/(np.amax(db) - np.amin(db)))
            elif cmcolor == 'tab20b':
                scattercolor = cm.tab20b(db)
    else:
        scattercolor = cm.Greys(np.ones(len(rep2d))*0.6)
    
    
    annotatemot = None
    annotate = None
    annotateloc = None
    if '--annotate' in sys.argv:
        # annotate motifs, for example
        annotateloc = 'inside'
        cmotiflist = sys.argv[sys.argv.index('--annotate')+1]
        annotateloc = sys.argv[sys.argv.index('--annotate')+2]
        minannosize = int(sys.argv[sys.argv.index('--annotate')+3])
        if cmotiflist == 'coloring':
            annotatemot = np.copy(coloringmot)
        elif cmotiflist == 'clustername':
            annotatemot = np.copy(clusters)
        if clustervisual != 'None':
            annotate = clusterlen >= minannosize
        else:
            annotate = np.ones(len(annotatemot)) == 1
    
   ### centresize: size of the dot is equal to number of rbps
    # centrecolor: all same size, only centre has color, rest is grey
    # add arrow for main factors
    # make multi-factor against each other plot


    ### define which proteins should get name text in 2d plot
    protext = []
    if '--locate_protein' in sys.argv:
        proteinnames = np.genfromtxt(sys.argv[sys.argv.index('--locate_protein')+1], dtype = str)
        to_mark_proteins = sys.argv[sys.argv.index('--locate_protein')+2].split(',')
        protext = []
        for tmp in to_mark_proteins:
            rpna = proteinnames[list(proteinnames[:,1]).index(tmp),3]
            rpnaloc = rep2d[list(prots).index(tmp)]
            protext.append([rpna, rpnaloc])
     
    
    colorfactor = []
    if '--color_factorplot' in sys.argv:
        numfactors = sys.argv[sys.argv.index('--color_factorplot')+1]
        pcafactor = sys.argv[sys.argv.index('--color_factorplot')+2] # whether to use pca or use original data from plotfeatures, 'None'
        combinedfactor = checkbool(sys.argv[sys.argv.index('--color_factorplot')+3]) # if true factors are colored in same umap plot
        normalize = checkbool(sys.argv[sys.argv.index('--color_factorplot')+4])
        if ',' in numfactors: 
            numfactors = np.array(numfactors.split(','), dtype = int)
        else:
            numfactors = np.arange(int(numfactors)+1, dtype = int)
        dfactors = reduce2d(plotfeatures, pcafactor, np.amax(numfactors)+1, None)
        
        colorfactor = dfactors[:, numfactors]
        colorfactor = colorfactor.T
        print sys.argv[sys.argv.index('--color_factorplot')+5], os.path.isfile(sys.argv[sys.argv.index('--color_factorplot')+5])
        if len(sys.argv) > sys.argv.index('--color_factorplot')+5 and os.path.isfile(sys.argv[sys.argv.index('--color_factorplot')+5]):
            proteinnames = np.genfromtxt(sys.argv[sys.argv.index('--color_factorplot')+5], dtype = str)
            protext = []
            for t in colorfactor:
                ams = np.argmin(t)
                amf = np.argmax(t)
                tmpmin = prots[ams]
                tmpmax = prots[amf]
                rpnamax = proteinnames[list(proteinnames[:,1]).index(tmpmax),3]
                rpnamin = proteinnames[list(proteinnames[:,1]).index(tmpmin),3]
                protext.append([[rep2d[ams], rpnamin], [rep2d[amf], rpnamax]])
        

    if '--figuresettings' in sys.argv:
        fwidth = float(sys.argv[sys.argv.index('--figuresettings')+1])
        fheight = float(sys.argv[sys.argv.index('--figuresettings')+2]) 
        falpha = float(sys.argv[sys.argv.index('--figuresettings')+3])
        fscatsize = float(sys.argv[sys.argv.index('--figuresettings')+4])
        faxisvis = checkbool(sys.argv[sys.argv.index('--figuresettings')+5])
        fconsize = float(sys.argv[sys.argv.index('--figuresettings')+6])
        fhead = float(sys.argv[sys.argv.index('--figuresettings')+7])
        anfontsize = int(sys.argv[sys.argv.index('--figuresettings')+8])
    else:
        fwidth = 10
        fheight = 10 
        falpha = 0.95
        fscatsize = 1.
        faxisvis = False
        fconsize = 3.  
        fhead = .5
        anfontsize = 10
    
    figsettings= [fwidth, fheight, falpha, fscatsize, faxisvis, fconsize, fhead, anfontsize]
    #scattercolor = ['grey' for i in range(len(scattercolor))]
    
    if '--savefig' in sys.argv:
        plotname = sys.argv[sys.argv.index('--savefig')+1]
    else:
        plotname = None
    
    if '--color_factorplot' in sys.argv:
        plotfactor(rep2d, colorfactor, numfactors, combinedfactor, protext, normalize, plotname, figsettings)
    
    if '--plotevo' in sys.argv:
        plotevo(rep2d, annotatemot, annotate, annotateloc, scattercolor, colorscheme, clustcentre, clustercomp, clustervisual, shape, connectmat, protext, plotname, figsettings)
    
    

    
    

    
