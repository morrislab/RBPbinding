#matrix_factorization.py 
#Identity_regression.py
import numpy as np
import scipy
import sys, os
import matplotlib
if "--savefig" in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from numpy.linalg import lstsq
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.linalg import orth
from numpy.linalg import inv
from scipy.linalg import svd
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
from scipy.spatial.distance import euclidean
from sklearn import decomposition
from sklearn.externals import joblib

if "--savematrix" in sys.argv:
    Wname = os.path.splitext(sys.argv[sys.argv.index("--savematrix")+1])[0]+".npz"
    print "Matrix outname is", Wname
    compute = False
elif "--recompute" in sys.argv:
    Wfactors = sys.argv[sys.argv.index("--recompute")+1]
    Wname = os.path.splitext(sys.argv[sys.argv.index("--recompute")+2])[0]+".npz"
    compute = True
    Facin = np.load(Wfactors)
    Htrain = Facin['factors']
    Wtrain = Facin['coefficients']
    seqfeatcompare = Facin['proteinfeatures']
    Zfeatcompare = Facin['Zfeatures']
    trainprots = Facin['trainprots']
    Ptrain = Facin['Ptrain']
    #HtrainP = Htrain[:, -len(Ptrain[0]):]
    #HtrainY = Htrain[:, :-len(Ptrain[0])]
    #print len(HtrainP[0]), len(HtrainY[0])
    #print np.sum(HtrainP, axis = 1), np.sum(HtrainY, axis = 1)
    #sys.exit()
else: 
    print "Please define whether to recompute results or save the generated matrix" 
    sys.exit()

if "--data" in sys.argv:
    profiles = sys.argv[sys.argv.index("--data")+1]
    proteinfeatures = sys.argv[sys.argv.index("--data")+2]
    print "Read in data..."
    print profiles
    #### get feature matrices and names
    Y = np.genfromtxt(profiles)[:,1:]
    Zfeatures = np.genfromtxt(profiles, dtype = str)[:,0]
    ynames = np.array(open(profiles, 'r').readline().strip().split()[1:])
    print proteinfeatures
    Y = np.nan_to_num(Y).T
    
    #determine name lines!
    if os.path.splitext(proteinfeatures)[1] == '.npz':
        pobj = np.load(proteinfeatures)
        P = pobj['features']
        sequencefeatures = pobj['kmers']
        protnames = pobj['protnames']
        pnames = pobj['expnames']
        pnames = np.append([pnames], [protnames], axis = 0)
    else:
        pobj = open(proteinfeatures, 'r').readlines()
        pnames = []
        sequencefeatures = []
        ph = 0
        for line in pobj:
            if line[0] == '#':
                ph +=1
                pnames.append(line.strip().split()[1:])
            else:
                sequencefeatures.append(line.split()[0])

        P = np.genfromtxt(proteinfeatures, skip_header = ph)[:,1:]
    
    P = P.T
    #sort after pnames the same order
    if len(np.intersect1d(ynames, pnames[0])) != 0:
        protnames = np.array(pnames[0])
    else:
        protnames = np.array(pnames[1])
    print protnames
    if not compute:
        indexes = []
        for i, ina in enumerate(protnames):
            if ina in ynames:
                indexes.append([i, list(ynames).index(ina)])
        indexes = np.array(indexes).T
    
    #print indexes, protnames[ :3], ynames[:3]
        Y = Y[indexes[1]]
        P = P[indexes[0]]
        protnames = protnames[indexes[0]]
        ynames = datnames = ynames[indexes[1]]
        mprot = np.shape(P)[0]
        nprobes = np.shape(Y)[1]
    else:
        indexes = []
        for i, ina in enumerate(ynames):
            if ina in trainprots:
                indexes.append(list(ynames).index(ina))
        Y = Y[indexes]
        ynames = ynames[indexes]
        
else:
    print "Please insert data or define simulation"
    sys.exit()

if compute:
    #comparesequence and Zscore features
    #Zfeatures, sequencefeatures, P, Y protnames, datnames
    print np.shape(Y)
    print np.shape(P)
    if not np.array_equal(seqfeatcompare, sequencefeatures):
        print 'Protein features differ', len(seqfeatcompare), len(sequencefeatures)
        norder = []
        opos = []
        print 'Num protfeatures', np.sum(P, axis = 1)
        Pnew = np.zeros(( len(seqfeatcompare), len(P)))
        sequencefeatures = list(sequencefeatures)
        for s, seqfeatc in enumerate(seqfeatcompare):
            if seqfeatc in sequencefeatures:
                norder.append(list(sequencefeatures).index(seqfeatc))
                opos.append(s)
        Pnew[opos] = P.T[norder]
        P = Pnew.T
        print 'Num protfeatures', np.sum(P, axis = 1)
        print len(protnames)
        
        sequencefeatures = seqfeatcompare
    if not np.array_equal(Zfeatcompare, Zfeatures):
        print 'Zscore features differ', len(Zfeatcompare), len(Zfeatures)
        norder = []
        for zfeatc in Zfeatcompare:
            norder.append(Zfeatures.index(zfeatc))
        Y = Y[:, norder]
        Zfeatures = Zfeatures[norder]
        
    
if "--cutP" in sys.argv:
    print "P contains binary values"
    P[P>0] = 1

if '--negY' in sys.argv:
    Y = -Y

if '--binaryY' in sys.argv:
    bcoff = float(sys.argv[sys.argv.index('--binaryY')+1])
    bctype = sys.argv[sys.argv.index('--binaryY')+2]
    if bctype == 'cut':
        if bcoff < 1:
            bcoff = int(np.shape(Y)[-1] * bcoff)
        else:
            bcoff = int(bcoff)
        for py, Yp in enumerate(Y):
            sort = np.argsort(-Yp)
            Y[py, sort[:bcoff]] = 1
            Y[py, sort[bcoff:]] = 0
    if bctype == 'sig':
        bcoffset = -np.sort(-Y.flatten())
        bcoff = bcoffset[int(len(bcoffset)*bcoff)]
        print 'Y lager than', bcoff
        for py, Yp in enumerate(Y):
            ymask = Y[py]>bcoff
            Y[py, ymask] = 1
            Y[py, ~ymask] = 0
        
if '--atleastkfeatures' in sys.argv:
    kmini = int(sys.argv[sys.argv.index('--atleastkfeatures')+1])
    sequencefeatures = np.array(sequencefeatures)[np.sum(P, axis = 0) > kmini]
    P = P[:,np.sum(P, axis = 0) > kmini]
    print 'All features removed that appeared less than', kmini, 'times!'

if '--minProteinfeatures' in sys.argv:
    Pmini = int(sys.argv[sys.argv.index('--atleastkfeatures')+1])
    mask = np.sum(P, axis = 1) >= pmini
    P = P[:,mask]
    Y = Y[mask]
    datnames = ynames = snames = protnames = protnames[mask]
    print 'All proteins removed that had less than', Pmini, 'features!'

if '--normP1' in sys.argv:
    print "normP1"
    P = (P.T/np.sum(np.absolute(P), axis = 1)).T

if '--cutY' in sys.argv:
    Y[Y<0.]= 0.

if '--preprocessingP' in sys.argv:
    P = (preprocessing.scale(P.T)).T
if '--preprocessingY' in sys.argv:
    Y = (preprocessing.scale(Y.T)).T

if '--minmaxY' in sys.argv:
    scaler = preprocessing.MinMaxScaler().fit(Y.T)
    Y = (scaler.transform(Y.T)).T
elif '--minmaxP' in sys.argv:
    scaler = preprocessing.MinMaxscaler().fit(P.T)
    P = (scaler.transform(P.T)).T
    
if '--normP2' in sys.argv:
    print "normP2"
    P = (P.T/np.sqrt(np.sum(P*P, axis = 1))).T

if '--normY2' in sys.argv:
    Y = (Y.T/np.sqrt(np.sum(Y*Y, axis = 1))).T

if '--normY1' in sys.argv:
    Y = (Y.T/np.sum(np.absolute(Y), axis = 1)).T

if "--quantile_norm" in sys.argv:
    print "Quantile normalization of Y, including test set"
    trainranks = np.argsort(Y, axis = 1).argsort(axis = 1)
    meanvalues = np.mean(np.sort(Y, axis =1), axis = 0)
    meanvalues = meanvalues/np.sqrt(np.sum(meanvalues*meanvalues))
    Y = meanvalues[trainranks]

if '--adjustP2Y' in sys.argv:
    Y = Y*float(len(P[-1]))/float(len(Y[-1]))

Y = np.nan_to_num(Y)
P = np.nan_to_num(P)
Y[np.isinf(Y)] = 0.
P[np.isinf(P)] = 0.
        
if "--trainingset" in sys.argv:
    numset = sys.argv[sys.argv.index("--trainingset")+1]
    trainfile = sys.argv[sys.argv.index("--trainingset")+2]
    print "trainset:", trainfile, numset
    tobj = open(trainfile, 'r')
    tlines = tobj.readlines()
    for i, tline in enumerate(tlines):
        if tline[:6]=="###Set" and tline.strip().split()[-1] == numset:
            trainprots = np.array([tlines[i+2].strip().split()])
            if tlines[i+3][:6] == '##Test':
                trainprots = np.append(trainprots, [tlines[i+2].strip().split()], axis = 0).T
                testprots = np.array([tlines[i+4].strip().split(), tlines[i+4].strip().split()]).T
            else:
                trainprots = np.append(trainprots, [tlines[i+3].strip().split()], axis = 0).T
                testprots = np.array([tlines[i+5].strip().split(), tlines[i+6].strip().split()]).T
            break
    ### sort into training and test set
    if len(np.intersect1d(trainprots[:,0], datnames)) !=0:
        trainprots = trainprots[:,0]
        testprots = testprots[:,0]
    else: 
        trainprots = trainprots[:,1]
        testprots = testprots[:,1]
    
    if compute:
        trainprots = trainprotsin[:]
        
    indexestrain = []
    indexestest = []
    for i, ina in enumerate(datnames):
        if ina in trainprots:
            indexestrain.append(i)
        if ina in testprots:
            indexestest.append(i)
    indexestrain = np.array(indexestrain)
    indexestest = np.array(indexestest)    
    
    Ytrain = Y[indexestrain]
    Ptrain = P[indexestrain]
    Ytest = Y[indexestest]
    Ptest = P[indexestest]
    trainprots = datnames[indexestrain]
    testprots = datnames[indexestest]
    print np.shape(Ytest), np.shape(Ptest)
    print np.shape(Ytrain), np.shape(Ptrain)

elif "--crossvalidation" in sys.argv:
    ncross = int(mprot/float(sys.argv[sys.argv.index("--crossvalidation")+1]))
    print sys.argv[sys.argv.index("--crossvalidation")+1]+"-fold crossvalidation" 
    Ptest = P[-ncross:]
    Ptrain = P[:mprot-ncross]
    Ytest = Y[-ncross:]
    Ytrain = Y[:mprot-ncross]
    if sim == False:
        trainprots = datnames[-ncross:]
        testprots = datnames[:mprot-ncross]
elif '--proteinlist' in sys.argv:
    plist = np.genfromtxt(sys.argv[sys.argv.index('--proteinlist')+1], dtype = str)
    pmask = np.isin(protnames, plist)
    if np.sum(pmask) == 0:
        print 'proteinlist does not match'
        print protnames, plist
        npnames = []
        for prna in protnames:
            if '|' in prna:
                prna = prna.split('|')[0]
            elif '_' in prna:
                prna = prna.split('_')[0]
            npnames.append(prna)
            pmask = np.isin(npnames, plist)
    print pmask
    protnames = protnames[pmask]
    datnames = np.copy(protnames)
    Ptest = P[pmask]
    print np.sum(Ptest, axis = 1)
    Ytrain = Ytest = Y[:]
    testprots = datnames
elif compute:
    datnames = np.copy(protnames)
    Ptest = P[:]
    print np.sum(Ptest, axis = 1)
    Ytrain = Ytest = Y[:]
    testprots = datnames
else:
    Ptrain = Ptest = P[:]
    Ytrain = Ytest = Y[:]
    trainprots = testprots = datnames

'''
if '--add_orthologs' in sys.argv:
    orthfeats = sys.argv[sys.argv.index("--add_orthologs")+1]
    print "Read in orthologs data..."
    Porth = np.load(features)['features']
    Orthnames = np.load(features)['protnames'] 
    orthchars = []
    
    for o, Orthname in enumerate(Orthnames):
        orthchars.append(Orthname.split('_', 1)[0])
    orthchars = np.array(orthchars)
    # add feature vectors of additional proteins from npz files and assign zscore to them 
    for t, trainprot in trainprots:
        if trainprot in orthchars:
            oloc = np.where(orthchars == trainprot)[0]
            Ptrain = np.append(Ptrain, Porth[oloc], axis = 0)
            Ytrain = np.append(Ytrain, [Ytrain[t] * len(oloc) ], axis = 0 ) 
'''
savefig = False
if "--savefig" in sys.argv:
    outname = sys.argv[sys.argv.index("--savefig")+1]
    print "Figurename", outname
    savefig = True
 
 
if compute == False: 
    Xtrain = np.append( Ytrain, Ptrain, axis = 1)

    import multiprocessing
    nj = multiprocessing.cpu_count()    
    if "--SVD" in sys.argv:
            perP = float(sys.argv[sys.argv.index("--SVD")+1])
            
            print "SVD\nproteinfeatures "+str(perP)
            Wtrain, Sigt, Ht = np.linalg.svd(Xtrain, full_matrices=False)
            
            dimP = min(len(Sigt), len(np.where(np.cumsum(Sigt)/np.sum(Sigt) < perP)[0])+1)
            print "new dim P", dimP

            Wtrain = Wtrain[:, :dimP]
            Sigt = Sigt[:dimP]
            Htrain = Ht[: dimP, :]
            
            print np.shape(Wtrain)
            print np.shape(Htrain)
            Wtrain = Wtrain*Sigt
    
    elif '--recommenderMF' in sys.argv:
            # may not work
            # IDEA: k-mer that correlates a lot with 7-mer has high entrance in covariance matrix
            # Now we can find latent space for k-mers and 7-mers. If k-mer important for 7-mer then correlation is high.
            from scipy import linalg
            from scipy.sparse.linalg import svds, eigs
            perP = float(sys.argv[sys.argv.index('--recommenderMF')+1])
            
            Xtrain = np.dot(Ptrain.T, Ytrain)
            print 'Recommend dimensions', np.shape(Xtrain)
            
            if perP > 1:
                dimP = int(perP)
                remethod = sys.argv[sys.argv.index('--recommenderMF')+2]
                if remethod == 'SVD':
                    print 'sparse SVD'
                    HtrainP, Sigt, HtrainY = svds(Xtrain, k = len(Ptrain))
                    print "new dim P", dimP
                elif remethod == 'NMF':
                    print 'NMF'
                    l1alpha = remethod = float(sys.argv[sys.argv.index('--recommenderMF')+3])
                    model = decomposition.NMF(n_components=dimP,  solver='mu', max_iter=5000, random_state=None, alpha=l1alpha, l1_ratio=1.)
                    model.fit(Xtrain)
                    HtrainY = model.components_
                    HtrainP = model.transform(Xtrain)
                    Sigt = np.ones(dimP) 
            else:
                print 'sparse SVD'
                HtrainP, Sigt, HtrainY = svds(Xtrain, k = len(Ptrain))
                dimP = min(len(Sigt), len(np.where(np.cumsum(Sigt)/np.sum(Sigt) < perP)[0])+1)
                print "new dim P", dimP

            HtrainP = HtrainP[:, :dimP]
            Sigt = Sigt[:dimP]
            HtrainY = HtrainY[: dimP, :]     
            
            HtrainP = HtrainP * np.sqrt(Sigt)
            HtrainY = (HtrainY.T * np.sqrt(Sigt)).T       
            
            lr = LinearRegression(fit_intercept = False).fit(HtrainP, Ptrain.T)
            Wtrain = lr.coef_
            
            #print linalg.norm(Ytrain - np.dot(Wtrain, HtrainY))
            #print linalg.norm(np.dot(Ytrain, Ytrain.T) - np.dot(np.dot(Wtrain, HtrainY), Ytrain.T))
            
            Htrain = np.append(HtrainY, HtrainP.T, axis = 1 )
            
            #print np.shape(Wtrain)
            #print np.shape(Htrain)

    elif '--recommenderEmbedding' in sys.argv:
            # maybe look if exponential family embeddings could be used!!! (External software)
            print 'not implemented'
            sys.exit()
            # needs to be implemented in pytorch: (StarSpace can potentially be used)
            # general idea: protein representation is given by the sum of k-mer latent-vectors
            # protein-latentvector multiplied by other protein latent-vectors are similarity of motifs
            # or protein latent-vector multiplied by 7-mer z-score latent vector is z-score for protein
            # Problem, learn latent entries for all protein k-mers and 7-mers? 
            
    
    elif '--regularizedSVD' in sys.argv:
            print 'not implemented'
            sys.exit()
            #https://github.com/yamitzky/Regularized-SVD
            #https://github.com/alabid/PySVD
            #https://github.com/metpallyv/MovieRecommendation

    elif '--graphregularizedNMF' in sys.argv:
            from scipy import linalg
            from joblib import Parallel, delayed
            import multiprocessing
            #from numba import jit
            num_cores = multiprocessing.cpu_count()
            # should be generalized and converted into Fortran code 
                    #k:          number of components
                    #W:          Adjacency matrix of the nearest neighboor matrix
                    #alpha:      regularization parameter
                    #init:       'random' initializes to a random W,H
                    #n_inits:    number of runs to make with different random inits (in order to avoid being stuck in local minima)
                    #n_jobs:     number of parallel jobs to run, when n_inits > 1
                    #tol:        stopping criteria
                    #max_iter:   stopping criteria


            def graphNMF(X, k=None, W=None, alpha=100, init='random', n_inits=1, tol=1e-3, max_iter=1000, update_rule = 'mu'):
                if k is None:
                    k = len(X[0])
                
                if W is None:
                    W = np.eye(len(X[0]))
                    D = np.eye(len(X[0]))
                    L = D - W
                else:
                    Wsum = np.sum(W, axis = 1)
                    D = np.eye(len(Wsum)) * Wsum
                    L = D - W
                
                if init == 'random':
                    U = np.random.random((X.shape[0], k))
                    V = np.random.random((X.shape[1], k))
                else:
                    print "Init not defined"
                    sys.exit()
                
                if update_rule == 'mu':
                    results = Parallel(n_jobs=num_cores)(delayed(gnmfcalc)(X, k, D, L, W, U, V, alpha, tol, max_iter) for i in range(n_inits))
                elif update_rule == 'cd':
                    results = Parallel(n_jobs=num_cores)(delayed(gnmfcalccd)(X, k, D, L, W, U, V, alpha, tol, max_iter, 0.0001) for i in range(n_inits))
                #results = gnmfcalc(X, k, D, L, U, V, alpha, init, n_inits, tol, max_iter)
                
                bresult = 1e10
                brun = None
                for r, result in enumerate(results):
                        if result[2] < bresult:
                            bresult = result[2]
                            brun = r
                U = results[brun][0]
                V = results[brun][1]
                return results[brun], np.linalg.norm(X - np.dot(U,V.T))
                
                
            def gnmfcalc(X, k, D, L, W, U, V, alpha, tol, max_iter):
                conv = False
                for x in range(max_iter):
                    #print x
                    XV = np.dot(X,V)
                    UVtV = np.dot(U, np.dot(V.T,V))

                    XtUpaWV = np.dot(X.T, U) + alpha * np.dot(W, V)
                    VUtUpaDV = np.dot(V, np.dot(U.T, U)) + alpha * np.dot(D, V)
                    
                    Un = np.multiply(U, np.nan_to_num(np.divide(XV, UVtV)))
                    Vn = np.multiply(V, np.nan_to_num(np.divide(XtUpaWV, VUtUpaDV)))
                    
                    e = linalg.norm(U-Un)
                    U, V = Un, Vn
                    #print e
                    if e < tol:
                        conv = True
                        break
                return U, V, e, conv
                
            #@jit(nopython = True, parallel = True)
            def gnmfcalccd(X, k, D, L, W, U, V, alpha, tol, max_iter, upalpha):
                tol = tol**2
                conv = False
                Uorig = np.copy(U)
                Vorig = np.copy(V)
                XV = np.dot(X,V)
                UVtV = np.dot(U, np.dot(V.T,V))
                Un = U + upalpha * (XV - UVtV)
                ean = np.sum((U - Un)**2) 
                for x in range(max_iter):
                    #print x, upalpha
                    XV = np.dot(X,V)
                    UVtV = np.dot(U, np.dot(V.T,V))

                    XtUpaWV = np.dot(X.T, U) + alpha * np.dot(W, V)
                    VUtUpaDV = np.dot(V, np.dot(U.T, U)) + alpha * np.dot(D, V)
                    
                    Un = U + upalpha * (XV - UVtV)
                    Vn = V + upalpha * (XtUpaWV - VUtUpaDV)
                    
                    e = np.sum((U - Un)**2)
                    U, V = Un, Vn
                    #print e
                    if e < tol:
                        conv = True
                        break
                    if e > ean:
                        upalpha = upalpha * 0.1
                        U = np.copy(Uorig)
                        V = np.copy(Vorig)
                        XV = np.dot(X,V)
                        UVtV = np.dot(U, np.dot(V.T,V))
                        Un = U + upalpha * (XV - UVtV)
                
                return None, None, e, conv
                
            print 'graphregularizedNMF'
            perP = float(sys.argv[sys.argv.index('--graphregularizedNMF')+1])
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)
            adjecmatfile = sys.argv[sys.argv.index('--graphregularizedNMF')+2]
            galpha = float(sys.argv[sys.argv.index('--graphregularizedNMF')+3])
            #n_inits= float(sys.argv[sys.argv.index('--graphregularizedNMF')+4])
            #max_iter = float(sys.argv[sys.argv.index('--graphregularizedNMF')+5])
            graphnames = open(adjecmatfile, 'r').readline().strip().split()[1:]
            adgraph = np.genfromtxt(adjecmatfile)
            if '--graphcut' in sys.argv:
                gcut = float(sys.argv[sys.argv.index('--graphcut')+1])
                adgraph[adgraph > gcut] = 1
                adgraph[adgraph <= gcut] = 0
                np.fill_diagonal(adgraph, 0)

            reorder = []
            for tp, trainp in enumerate(trainprots):
                reorder.append(graphnames.index(trainp))
            adgraph = adgraph[reorder]
            adgraph = adgraph[:, reorder]
            
            print np.shape(Xtrain.T)
            model, erro = graphNMF(Xtrain.T, k=dimP, W=adgraph, alpha=galpha, init='random', n_inits=10, max_iter=5000, update_rule = 'cd')
            print 'Reconstruction error', erro
            Htrain = model[0].T
            Wtrain = model[1]
            
    elif "--NMF" in sys.argv:
            # use NMF when less features than data?
            print "NMF..."
            perP = float(sys.argv[sys.argv.index('--NMF')+1])
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)
            if '--nmfalpha' in sys.argv:
                l1alpha = float(sys.argv[sys.argv.index('--nmfalpha')+1])
                l1ratio = float(sys.argv[sys.argv.index('--nmfalpha')+2])
            if '--changesolver' in sys.argv:
                solver = 'mu'
            else:
                solver = 'cd'
            if '--changeloss' in sys.argv:
                beta_loss='kullback-leibler'
            else:
                beta_loss='frobenius'
            model = decomposition.NMF(n_components=dimP,  solver=solver, max_iter=5000, random_state=None, alpha=l1alpha, l1_ratio=l1ratio)
            model.fit(Xtrain)
            Htrain = model.components_
            Wtrain = model.transform(Xtrain)

    elif '--regularizedPCA' in sys.argv:
            print ' Regularized PcA'
            perP = float(sys.argv[sys.argv.index('--regularizedPCA')+1])
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)  
            
            l1alpha = float(sys.argv[sys.argv.index('--regularizedPCA')+2])
            l2alpha = float(sys.argv[sys.argv.index('--regularizedPCA')+3])
            
            # regularized PCA, ridge_alpha is L2 norm. alpha is L1 norm
            model = decomposition.SparsePCA(n_components=dimP, alpha=l1alpha, ridge_alpha=l2alpha, max_iter=10000, tol=1e-08, method='lars', n_jobs=nj, U_init=None, V_init=None, verbose=False, random_state=None)
            model.fit(Xtrain)
            
            
            Htrain = model.components_
            Wtrain = model.transform(Xtrain)

    elif '--sparseSVD' in sys.argv:
            from scipy.sparse.linalg import svds, eigs
            
            perP = float(sys.argv[sys.argv.index('--sparseSVD')+1])
            Wtrain, Sigt, Ht = svds(Xtrain, k=len(Ptrain))
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)
            
            dimP = min(len(Sigt), len(np.where(np.cumsum(Sigt)/np.sum(Sigt) < perP)[0])+1)
            
            Wtrain = Wtrain[:, :dimP]
            Sigt = Sigt[:dimP]
            Htrain = Ht[: dimP, :]
            
            print np.shape(Wtrain)
            print np.shape(Htrain)
            Wtrain = Wtrain*Sigt

    elif '--DictionaryLearning' in sys.argv:
            print 'dictionary learning'
            #(U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || U ||_1
            #    (U,V)
            #   with || V_k ||_2 = 1 for all  0 <= k < n_components
            perP = float(sys.argv[sys.argv.index('--DictionaryLearning')+1])
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)
            model = decomposition.DictionaryLearning(n_components=dimP, alpha=1., max_iter=10000, tol=1e-08, fit_algorithm='lars', transform_algorithm='lasso_lars', transform_n_nonzero_coefs=None, transform_alpha=None, n_jobs=nj, code_init=None, dict_init=None, verbose=False, split_sign=False, random_state=None)
            model.fit(Xtrain)
            Htrain = model.components_
            Wtrain = model.transform(Xtrain)  
            
    elif '--FactorAnalysis' in sys.argv:
            print 'FactorAnalysis, a classical statistical model. The matrix W is sometimes called the "factor loading matrix"'
            #Both models essentially estimate a Gaussian with a low-rank covariance matrix. Because both models are probabilistic they can be integrated in more complex models, e.g. Mixture of Factor Analysers. One gets very different models (e.g. FastICA) if non-Gaussian priors on the latent variables are assumed
            #The main advantage for Factor Analysis (over PCA is that it can model the variance in every direction of the input space independently (heteroscedastic noise)
            #The initial guess of the noise variance for each feature. If None, it defaults to np.ones(n_features)
            perP = float(sys.argv[sys.argv.index('--FactorAnalysis')+1])
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)
            model = decomposition.FactorAnalysis(n_components=dimP, tol=0.01, copy=True, max_iter=1000, noise_variance_init=None, svd_method='randomized', iterated_power=3, random_state=0)
            model.fit(Xtrain)
            Htrain = model.components_
            Wtrain = model.transform(Xtrain)
            print 'Noise variance', model.noise_variance_

    elif '--FastIcA' in sys.argv:
            # could be used for Y
            print 'Independent component analysis (ICA) is used to estimate sources given noisy measurements. Imagine 3 instruments playing simultaneously and 3 microphones recording the mixed signals. ICA is used to recover the sources ie. what is played by each instrument. Importantly, PCA fails at recovering our instruments since the related signals reflect non-Gaussian processes.'
            #### Define your own function f.e.: def g(x):
            perP = float(sys.argv[sys.argv.index('--FastIcA')+1])                                    #return x**3, 3x**2 
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)
            model = decomposition.FastICA(n_components=dimP, algorithm='parallel', whiten=True, fun='logcosh', fun_args=None, max_iter=10000, tol=0.0001, w_init=None, random_state=None)
            model.fit(Xtrain)
            Htrain = model.components_
            Wtrain = model.transform(Xtrain)

    elif '--LDA' in sys.argv:
            #LDA could also be more valuable for Ys
            print 'Latent dirichlet allocation'
            perP = float(sys.argv[sys.argv.index('--LDA')+1])
            if perP > 1:
                dimP = int(perP)
            else:
                dimP = int(len(Ptrain)*perP)
            model = decomposition.LatentDirichletAllocation(n_components=dimP, doc_topic_prior=None, topic_word_prior=None, learning_method=None, learning_decay=0.7, learning_offset=10.0, max_iter=10, batch_size=128, evaluate_every=-1, total_samples=1000000.0, perp_tol=0.1, mean_change_tol=0.001, max_doc_update_iter=100, n_jobs=1, verbose=0, random_state=None, n_topics=None)
            model.fit(Xtrain)
            Htrain = model.components_
            Wtrain = model.transform(Xtrain)
    elif '--customMF' in sys.argv:
        custom = 1
        # insert lists of k-mers that can be combined
                # insert lists of l-mers that can be combined
                # number of p-feature and Z-score groups
                # objective functions, look into PCA and generalization of MF methods
                # constraints for factor values

    else:
            print "Please define method to reduce dimensionality"
            sys.exit()

    if '--savefactors' in sys.argv:
        print 'Factors saved'
        np.savez_compressed(os.path.splitext(Wname)[0]+'-factors.npz', factors = Htrain, coefficients = Wtrain, proteinfeatures = sequencefeatures, Zfeatures = Zfeatures, trainprots = trainprots, Ptrain = Ptrain)
 
# Evaluation of prediction
if '--lsq' in sys.argv:
    def transform_reconstruct(Invec, Infac, RecFac):
        lr = LinearRegression(fit_intercept = False).fit(Infac, Invec)
        Wcoef = lr.coef_
        VectRect = np.dot(Wcoef, RecFac)
        return Wcoef, VectRect
    
elif '--nnlsq' in sys.argv:
    def transform_reconstruct(Invec, Infac, RecFac):
        Wcoef = []
        for p, ptst in enumerate(Invec.T):
            Wcoef.append(nnls(Infac, ptst)[0])
        Wcoef = np.array(Wcoef)
        VectRect = np.dot(Wcoef, RecFac)
        return Wcoef, VectRect

elif '--L1reg' in sys.argv:
    alph = float(sys.argv[sys.argv.index('--L1reg')+1])
    pos = bool(sys.argv[sys.argv.index('--L1reg')+2])
    def transform_reconstruct(Invec, Infac, RecFac):
        lr = Lasso(alpha = alph, fit_intercept = False, positive = pos).fit(Infac, Invec)
        Wcoef = lr.coef_
        VectRect = np.dot(Wcoef, RecFac)
        return Wcoef, VectRect

HtrainP = Htrain[:, -len(Ptrain[0]):]
HtrainY = Htrain[:, :-len(Ptrain[0])]
Wtest, Ypred = transform_reconstruct(Ptest.T, HtrainP.T, HtrainY)


#elif '--reverseICA' in sys.argv:
    ##### Need another method for ICA: Could be that algorithm gets stuck in local minimum. Doesn't look like ICA learns something useful here. 
    #### ICA needs another method to fit the coefficients of the test protein. Look into the distribution of coefficients for the training proteins. 
    #### ICA generates coefficient profiles with only one, 1 or -1 and the rest 0.
    ##print np.shape(Wtrain)
    ##figw = plt.figure(1)
    ##figv= plt.figure(2)
    ##wl = len(Wtrain[0])
    ##vl = len(Wtrain)
    ##for i, w in enumerate(Wtrain):
        ##axw = figw.add_subplot(1,11,i+1)
        ##axw.hist(w, bins=int(float(wl)/3.), alpha = 0.5+i*0.05)
        ##if i == 10:
            ##break
    ##for i, v in enumerate(Wtrain.T):
        ##axv = figv.add_subplot(1,11,i+1)
        ##axv.hist(v, bins=int(float(vl)/3.), alpha = 0.5+i*0.05)
        ##if i == 10:
            ##break
    ##plt.show()
    #sys.exit()

#elif '--customcoeff' in sys.argv:
    #def sparsecoefficients(maxcoeff, stepsize, seqfeatures, seqfactors, objectivem, positcoeff = False):
        ## define objective function to maximize
        #if objectivem == 'pearson':
            #def obf(pprofile, predprofile):
                #return 1. - np.nan_to_num(pearsonr(pprofile, predprofile)[0])
        #elif objectivem == 'squared_error':
            #def obf(pprofile, predprofile):
                #return np.sqrt(np.mean((pprofile - predprofile)**2))
        #elif objectivem == 'absolute_error': 
            #def obf(pprofile, predprofile):
                #return np.mean(np.absolute(pprofile - predprofile))
        #elif objectivem == 'truepositive':
            #def obf(pprofile, predprofile):
                #poskmer = predprofile > np.amin(predprofile)
                #return np.sqrt(np.mean((pprofile[poskmer] - predprofile[poskmer])**2))
        #elif objectivem == 'trueposiveboth':
            #def obf(pprofile, predprofile):
                #poskmer = predprofile > np.amin(predprofile)
                #poskmer = poskmer * posfactor
                #return np.sqrt(np.mean((pprofile[poskmer] - predprofile[poskmer])**2))            
        #elif 'expected-realdistance' :
            #def obf(pprofile, predprofile):
                #return
        #elif objectivem == 'mutinfo': 
            #def obf(pprofile, predprofile):
                ## do histogram for each. 
                ## look up indiviual proability
                #return 1. - np.mean(np.absolute(pprofile - predprofile))
        
        #pcoeff = np.zeros(len(seqfactors))
        #oldbest = obf(np.zeros(len(seqfeatures)), seqfeatures)
        #print oldbest
        #while np.sum(pcoeff) <= maxcoeff:
            #if positcoeff:
                #testcoeff = np.ones((len(pcoeff), len(pcoeff))) * pcoeff + np.eye(len(pcoeff)) * stepsize
                #temppre = np.dot(testcoeff.T, seqfactors)
                #teststat = np.zeros(len(pcoeff))
                #for t, testc in enumerate(temppre):
                    #teststat[t] = obf(testc, seqfeatures)
            #else:
                #testcoeff = np.ones((len(pcoeff), len(pcoeff))) * pcoeff + np.eye(len(pcoeff)) * stepsize
                #testcoeff = np.append(testcoeff, np.ones((len(pcoeff), len(pcoeff))) * pcoeff - np.eye(len(pcoeff)) * stepsize, axis = 1)
                #temppre = np.dot(testcoeff.T, seqfactors)
                #print np.shape(temppre)
                #teststat = np.zeros(len(pcoeff)*2)
                #for t, testc in enumerate(temppre):
                    #teststat[t] = obf(testc, seqfeatures)
            #newbest = np.argsort(teststat)[0]
            #print newbest, teststat[newbest]
            #if teststat[newbest] <= oldbest:
                #pcoeff = testcoeff[:, newbest]
                #oldbest = teststat[newbest]
                #print oldbest, np.nonzero(pcoeff)[0]
            #else:
                #break
        #print pcoeff    
        #return pcoeff
    #HtrainP = Htrain[:, -len(Ptrain[0]):]
    #HtrainY = Htrain[:, :-len(Ptrain[0])]
    #Wtest = []
    #for p, ptst in enumerate(Ptest):
        #Wtest.append(sparsecoefficients(4, 1, ptst, HtrainP, 'truepositiveboth'))
        ##Wtest.append(sparsecoefficients(4, .2, Ytest[p], HtrainY, 'squared_error'))
    #Wtest = np.array(Wtest)
    #Ypred = np.dot(Wtest, HtrainY)
    ## maximum as many coefficients as domains + 1
    ## minimum of certain value for coefficients. 

    #print reconstlist
    
if '--lsq' in sys.argv or '--nnlsq' in sys.argv or '--L1reg' in sys.argv or '--customcoeff' in sys.argv:    
    if '--backtransform' in sys.argv:
        Btrans = sys.argv[sys.argv.index("--backtransform")+1]
        if os.path.splitext(Btrans)[1] == ".npz":
            Tp = np.load(Btrans)["Tp"]
        elif os.path.splitext(Btrans)[1] == ".txt":
            Tp = np.genfromtxt(Btrans, astype = str)
            embedfeats = Tp[:,0]
            Tp = Tp[:,1:].T.atype(float)
        elif os.path.splitext(Btrans)[1] == ".pkl":
            Tp = joblib.load(Btrans)[-1] 
        Ypred = np.dot(Ypred,Tp)
        Trainzscores = sys.argv[sys.argv.index("--backtransform")+2]
        Trainprofiles = np.genfromtxt(Trainzscores)[:,1:]
        Tnames = np.array(open(Trainzscores, 'r').readline().strip().split()[1:])
        Trainprofiles = np.nan_to_num(Trainprofiles).T
        Ytest = []
        for tprot in testprots:
            Ytest.append(Trainprofiles[Tnames == tprot][0])
        Ytrain = []
        for tprot in trainprots:
            Ytrain.append(Trainprofiles[Tnames == tprot][0])
        Ytest = np.array(Ytest)
        Ytrain = np.array(Ytrain)


        
        
    
    if '--save_reconstruction_error' in sys.argv:
        recmeasure = sys.argv[sys.argv.index('--save_reconstruction_error')+1]
        predicted_protein = np.dot( Wtest, HtrainP)
        def reconsterror(rmeasure, pprofile, predprofile):
            recorrectness = []
            #print '\n', len(np.nonzero(pprofile > 0)[0])
            #print len(np.nonzero(predprofile > 0)[0])
            #print len(np.nonzero((predprofile > 0)*(pprofile > 0))[0])
            if rmeasure == 'all':
                repearson = pearsonr(pprofile, predprofile)[0]
                reprecision = float(np.sum((pprofile > 0)*(predprofile > 0)))/max(1.,float(np.sum((predprofile > 0))))
                recthresh = np.sort(predprofile)[-int(0.1*float(len(predprofile)))]
                rectop10 = float(np.sum((pprofile > 0)*(predprofile > recthresh)))/max(1.,float(np.sum((predprofile > recthresh))))
                reverseprofile = np.ones(len(pprofile))- (pprofile > 0)
                reversepredictioin = np.ones(len(predprofile))- (predprofile > 0)
                renegrate = np.sum(reverseprofile*reversepredictioin)/np.sum(reversepredictioin)
                resqerror = 1. - np.sum((pprofile - predprofile)**2)
                recorrectness = [repearson, reprecision, rectop10, renegrate, resqerror]
            elif rmeasure == 'pearson':
                recorrectness.append(pearsonr(pprofile, predprofile)[0])
            elif rmeasure == 'precision':
                recorrectness.append(float(np.sum((pprofile > 0)*(predprofile > 0)))/float(np.sum((predprofile > 0))))
            elif rmeasure == 'precisiontop10' in sys.argv:
                recthresh = np.sort(predprofile)[-int(0.1*float(len(predprofile)))]
                recorrectness.append(float(np.sum((pprofile > 0)*(predprofile > recthresh)))/float(np.sum((predprofile > recthresh))))
            elif rmeasure == 'negativerate':
                reverseprofile = np.ones(len(pprofile))- (pprofile > 0)
                reversepredictioin = np.ones(len(predprofile))- (predprofile > 0)
                recorrectness.append(np.sum(reverseprofile*reversepredictioin)/np.sum(reversepredictioin))
            elif rmeasure == 'squared_error':
                recorrectness.append(1. - np.sum((pprofile - predprofile)**2))
            return recorrectness
        reconstlist = []
        for i in range(len(testprots)):
            reconstlist.append(reconsterror(recmeasure, predicted_protein[i], Ptest[i]))
        reconstlist = np.array(reconstlist)
        np.savetxt(os.path.splitext(Wname)[0]+"-reconsterror_"+recmeasure+".dat", np.append(testprots.reshape((len(testprots),1)), reconstlist.astype(str), axis=1), fmt = '%s')
        
    
    
    if '--savetestcorrelation' in sys.argv:
        pearsont = np.zeros((len(testprots),2))
        for i in range(len(testprots)):
            pearsont[i,0] = pearsonr(Ypred[i], Ytest[i])[0]
            pearsont[i,1] = spearmanr(Ypred[i], Ytest[i])[0]
            print testprots[i], pearsont[i,0], pearsont[i,1]
        np.savetxt(os.path.splitext(Wname)[0]+"-testset_profile_pcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')
        
        if '--run_analysis' in sys.argv:
            def featureoverlap(ftest, ftrain):
                testmin = np.amin(ftest)
                trainmin = np.amin(ftrain)
                ftest = ftest > testmin
                ftrain = ftrain > trainmin
                olapfeat = np.sum(ftest*ftrain)
                return olapfeat
            
            anafile = open(os.path.splitext(Wname)[0]+'_analysis.dat', 'w')
            for i in range(len(testprots)):
                if pearsont[i,0] < 0.6:
                    hastwin = 0
                    featlap = []
                    samelap = []
                    for j in range(len(trainprots)):
                        ptotrain = pearsonr(Ypred[i], Ytrain[j])[0] 
                        ptwtrain = pearsonr(Ytest[i], Ytrain[j])[0]
                        if ptotrain > 0.55:
                            samelap.append([trainprots[j], featureoverlap(Ptest[i], Ptrain[j]), ptotrain])
                        if ptwtrain > 0.55:
                            hastwin += 1
                            featlap.append([trainprots[j],featureoverlap(Ptest[i], Ptrain[j]), ptwtrain])
                        
                        
                    anafile.write('\nFalse: '+testprots[i]+' '+str( pearsont[i,0])+'\n')
                    anafile.write("Protein features "+str( featureoverlap(Ptest[i], Ptest[i]))+'\n')
                    if hastwin == 0:
                        anafile.write('# No similar protein!\n')
                    else:
                        anafile.write('# Similar proteins in trainingset\n')
                        for dto in featlap:
                            anafile.write(dto[0]+' '+str(dto[1])+' '+str(dto[2])+'\n')
                    anafile.write("# Predicted Similarity\n")
                    for dto in samelap:
                        anafile.write(dto[0]+' '+str(dto[1])+' '+str(dto[2])+'\n')
        if '--reconstruction_error' in sys.argv and '--analysereconstruction' in sys.argv:
            for i in range(len(testprots)):
                print testprots[i], pearsont[i,0], ' '.join(reconstlist[i].astype(str))
            print 'Total correlations'
            for rmes in reconstlist.T:
                print pearsonr(pearsont[:,0], rmes)[0]
                plt.scatter(pearsont[:,0], rmes)
            plt.show()
        
    if '--savetopintersection' in sys.argv:
        pearsont = np.zeros((len(testprots),2))
        perctop = int(len(Ypred[0])*0.01)
        for i in range(len(testprots)):
            pearsont[i,0] = len(np.intersect1d(np.argsort(Ypred[i])[-perctop:], np.argsort(Ytest[i])[-perctop:]))/float(perctop)
            pearsont[i,1] = len(np.intersect1d(np.argsort(Ypred[i])[-100:], np.argsort(Ytest[i])[-100:]))/(100.)       
            print testprots[i], pearsont[i,0], pearsont[i,1]
        np.savetxt(os.path.splitext(Wname)[0]+"-testset_profile_topintersection.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')
    
    if '--saveclosestcorrelation' in sys.argv:
        pearsont = []
        for i in range(len(testprots)):
            pearsi = []
            for j in range(len(trainprots)):
                pearsi.append(pearsonr(Ypred[i], Ytrain[j])[0])
            pearsont.append([trainprots[np.argmax(pearsi)], pearsi[np.argmax(pearsi)]] )
        #print pearsont
        np.savetxt(os.path.splitext(Wname)[0]+"-directfit_highcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), np.array(pearsont).astype(str), axis=1), fmt = '%s')
    
    if "--savetestprofiles" in sys.argv:
            print "Test Z Profiles saved"
            np.savez_compressed(os.path.splitext(Wname)[0]+"-testset_profiles.npz", names=testprots, profiles=Ypred.T, kmers = Zfeatures)
    
    if "--savetestprofileslight" in sys.argv:
           tk = []
           for t, tprot in enumerate(testprots):
               ksort = np.argsort(-Ypred[t])[:100]
               tk.append(Zfeatures[ksort])
           np.savetxt(os.path.splitext(Wname)[0]+"-testset_profiles_light.txt", np.array(tk).T, fmt= '%s', header=' '.join(testprots))

    if '--save_reconstructionP' in sys.argv:
        WtestZ, YtestZtoP = transform_reconstruct(Ytest.T, HtrainY.T, HtrainP)
        np.savez_compressed(os.path.splitext(Wname)[0]+"-testset_PreconstructionZ.npz", names=testprots, profiles=YtestZtoP, kmers = sequencefeatures)
        
    if '--save_reconstructionPtoP' in sys.argv:
        YtestPtoP = np.dot(Wtest, HtrainP)
        np.savez_compressed(os.path.splitext(Wname)[0]+"-testset_PreconstructionP.npz", names=testprots, profiles= YtestPtoP, kmers = sequencefeatures)
    
    if '--save_reconstructionPZtoP' in sys.argv:
        YtrainPZtoP = np.dot(Wtrain, HtrainP) # Wtrain equal WtestPZ in case train and test set equal
        np.savez_compressed(os.path.splitext(Wname)[0]+"-testset_PreconstructionPZ.npz", names=testprots, profiles = YtrainPZtoP, kmers = sequencefeatures)
        
    if '--savetopcorrelation' in sys.argv:
                pearsont = np.zeros((len(testprots),2))
                perctop = int(len(Ypred[0])*0.01)
                for i in range(len(testprots)):
                    weight100 = np.zeros(len(Ypred[i]))
                    weight1p = np.zeros(len(Ypred[i]))
                    w1 = np.argsort(Ypred[i])[-perctop:]
                    w2 = np.argsort(Ytest[i])[-perctop:]
                    weight1p[w1] += 0.5
                    weight1p[w2] += 0.5
                    weight100[w1[-100:]] +=0.5
                    weight100[w2[-100:]] +=0.5
                    mask1p = weight1p > 0
                    mask100 = weight100 > 0
                    weight1p = weight1p[mask1p]
                    weight100 = weight100[mask100]
                    pearsont[i,0] = np.sum((Ypred[i][mask1p]-np.mean(Ypred[i][mask1p]))*weight1p*(Ytest[i][mask1p]-np.mean(Ytest[i][mask1p])))/(np.sum(weight1p)*np.std(Ypred[i][mask1p])*np.std(Ytest[i][mask1p]))
                    pearsont[i,1] = np.sum((Ypred[i][mask100]-np.mean(Ypred[i][mask100]))*weight100*(Ytest[i][mask100]-np.mean(Ytest[i][mask100])))/(np.sum(weight100)*np.std(Ypred[i][mask100])*np.std(Ytest[i][mask100]))
                np.savetxt(os.path.splitext(Wname)[0]+"-testset_profile_topcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')    

    if '--savelatenttest' in sys.argv:
        np.savez_compressed(os.path.splitext(Wname)[0]+"-testset_Platent.npz", names=testprots, profiles=Wtest)

if '--predS' in sys.argv:
    smeasure = sys.argv[sys.argv.index('--predS')+1]
    print np.shape(Wtrain), np.shape(Wtest)
    
    if smeasure == 'dot': 
        Spred = np.dot(Ypred, Ytrain.T)

    elif smeasure == 'pearson':
        Ytrain = (preprocessing.scale(Ytrain.T)).T
        Ypred = preprocessing.scale(Ypred)
        Spred = np.zeros((len(Wtest),len(Wtrain)))
        for da in range(len(Ypred)):
            for dq in range(len(Ytrain)):
                Spred[da,dq] = pearsonr(Ypred[da], Ytrain[dq])[0]
        #Spred = np.dot(Ypred, Ytrain.T)

    elif smeasure == 'euclidean':
        Ytrain = (preprocessing.scale(Ytrain.T)).T
        Ypred = preprocessing.scale(Ypred)
        Spred = np.zeros((len(Ypred),len(Ytrain)))
        for da in range(len(Ypred)):
            for dq in range(len(Ytrain)):
                Spred[da,dq] = euclidean(Ypred[da], Ytrain[dq])

    elif smeasure == 'pearsonreconstructed':

        WtrainZ, YtrainZ = transform_reconstruct(Ytrain.T, HtrainY.T, HtrainY)
        YtrainZ = (preprocessing.scale(YtrainZ.T)).T
        Ypred = preprocessing.scale(Ypred)
        Spred = np.zeros((len(Wtest),len(Wtrain)))
        for da in range(len(Ypred)):
            for dq in range(len(YtrainZ)):
                Spred[da,dq] = pearsonr(Ypred[da], YtrainZ[dq])[0]
        
    elif smeasure == 'euclideanreconstructed':
    
        WtrainZ, YtrainZ = transform_reconstruct(Ytrain.T, HtrainY.T, HtrainY)
        YtrainZ = (preprocessing.scale(YtrainZ.T)).T
        Ypred = preprocessing.scale(Ypred)
        Spred = np.zeros((len(Wtest),len(Wtrain)))
        for da in range(len(Ypred)):
            for dq in range(len(YtrainZ)):
                Spred[da,dq] = euclidean(Ypred[da], YtrainZ[dq])
                
    elif smeasure == 'pearsonPtoZreconstructed':

        WtrainP, YtrainPtoZ = transform_reconstruct(Ptrain.T, HtrainP.T, HtrainY)
        YtrainPtoZ = (preprocessing.scale(YtrainPtoZ.T)).T
        Ypred = preprocessing.scale(Ypred)
        Spred = np.zeros((len(Wtest),len(Wtrain)))
        for da in range(len(Ypred)):
            for dq in range(len(YtrainPtoZ)):
                Spred[da,dq] = pearsonr(Ypred[da], YtrainPtoZ[dq])[0]
        
    elif smeasure == 'euclideanPtoZreconstructed':
    
        WtrainP, YtrainPtoZ = transform_reconstruct(Ptrain.T, HtrainP.T, HtrainY)
        YtrainPtoZ = (preprocessing.scale(YtrainPtoZ.T)).T
        Ypred = preprocessing.scale(Ypred)
        Spred = np.zeros((len(Wtest),len(Wtrain)))
        for da in range(len(Ypred)):
            for dq in range(len(YtrainPtoZ)):
                Spred[da,dq] = euclidean(Ypred[da], YtrainPtoZ[dq])
    
    elif smeasure == 'dotlatent': 
        Spred = np.dot(Wtest, Wtrain.T)
    
    elif smeasure == 'pearsonlatent':
        Wtrain = (preprocessing.scale(Wtrain.T)).T
        Wtest = (preprocessing.scale(Wtest.T)).T
        Spred = np.zeros((len(Wtest),len(Wtrain)))
        for da in range(len(Wtest)):
            for dq in range(len(Wtrain)):
                Spred[da,dq] = pearsonr(Wtest[da], Wtrain[dq])[0]
    
    elif smeasure == 'euclideanlatent':
        Wtrain = (preprocessing.scale(Wtrain.T)).T
        Wtest = (preprocessing.scale(Wtest.T)).T
        Spred = np.zeros((len(Wtest),len(Wtrain)))
        for da in range(len(Wtest)):
            for dq in range(len(Wtrain)):
                Spred[da,dq] = euclidean(Wtest[da], Wtrain[dq])
        scaler = preprocessing.MinMaxScaler().fit(Spred.T)
        Spred = 1. - (scaler.transform(Spred.T)).T
    
    elif smeasure == 'pearsonlatentZ':

        WtrainZ, YtrainZ = transform_reconstruct(Ytrain.T, HtrainY.T, HtrainY)
        WtrainZ = (preprocessing.scale(WtrainZ.T)).T
        Wtest = (preprocessing.scale(Wtest.T)).T
        Spred = np.zeros((len(Wtest),len(WtrainZ)))
        for da in range(len(Wtest)):
            for dq in range(len(WtrainZ)):
                Spred[da,dq] = pearsonr(Wtest[da], WtrainZ[dq])[0]
        
    elif smeasure == 'euclideanlatentZ':
    
        WtrainZ, YtrainZ = transform_reconstruct(Ytrain.T, HtrainY.T, HtrainY)
        WtrainZ = (preprocessing.scale(WtrainZ.T)).T
        Wtest = (preprocessing.scale(Wtest.T)).T
        Spred = np.zeros((len(Wtest),len(WtrainZ)))
        for da in range(len(Wtest)):
            for dq in range(len(WtrainZ)):
                Spred[da,dq] = euclidean(Wtest[da], WtrainZ[dq])
                
    elif smeasure == 'pearsonlatentP':

        WtrainP, YtrainPtoZ = transform_reconstruct(Ptrain.T, HtrainP.T, HtrainY)
        WtrainP = (preprocessing.scale(WtrainP.T)).T
        Wtest = (preprocessing.scale(Wtest.T)).T
        Spred = np.zeros((len(Wtest),len(WtrainP)))
        for da in range(len(Wtest)):
            for dq in range(len(WtrainP)):
                Spred[da,dq] = pearsonr(Wtest[da], WtrainP[dq])[0]
        
    elif smeasure == 'euclideanlatentP':
    
        WtrainP, YtrainPtoZ = transform_reconstruct(Ptrain.T, HtrainP.T, HtrainY)
        WtrainP = (preprocessing.scale(WtrainP.T)).T
        Wtest = (preprocessing.scale(Wtest.T)).T
        Spred = np.zeros((len(Wtest),len(WtrainP)))
        for da in range(len(Wtest)):
            for dq in range(len(WtrainP)):
                Spred[da,dq] = euclidean(Wtest[da], WtrainP[dq])  
    
    
    
    elif smeasure == 'pearsonlatentPtoZ':
        WtrainP, YtrainPtoZ = transform_reconstruct(Ptrain.T, HtrainP.T, HtrainY)
        WtrainZ, YtrainZ = transform_reconstruct(Ytrain.T, HtrainY.T, HtrainY)
        WtrainP = (preprocessing.scale(WtrainP.T)).T
        Wtest = (preprocessing.scale(Wtest.T)).T
        WtrainZ = (preprocessing.scale(WtrainZ.T)).T
        SpredZ = np.zeros((len(WtrainZ),len(WtrainZ)))
        for da in range(len(WtrainZ)):
            for dq in range(len(WtrainZ)):
                SpredZ[da,dq] = pearsonr(WtrainP[da], WtrainP[dq])[0]
        SpredP = np.zeros((len(Wtest),len(WtrainZ)))
        for da in range(len(Wtest)):
            for dq in range(len(WtrainP)):
                SpredP[da,dq] = pearsonr(Wtest[da], WtrainP[dq])[0]
        Spred = np.zeros((len(SpredP),len(SpredZ)))
        for da in range(len(SpredP)):
            for dq in range(len(SpredZ)):
                SpredP[da,dq] = pearsonr(SpredP[da], SpredZ[dq])[0]
        
        
    
    #elif smeasure == 'pearsonlatentPtoZtestplot':
        #WtrainP, YtrainPtoZ = transform_reconstruct(Ptrain.T, HtrainP.T, HtrainY)
        #WtrainZ, YtrainZ = transform_reconstruct(Ytrain.T, HtrainY.T, HtrainY)
        #WtrainP = (preprocessing.scale(WtrainP.T)).T
        #WtestP = (preprocessing.scale(Wtest.T)).T
        #WtestZ, YtestZ = transform_reconstruct(Ytest.T, HtrainY.T, HtrainY)
        #WtestZ = (preprocessing.scale(WtestZ.T)).T
        #SpredP = np.zeros((len(WtestP),len(WtrainP)))
        ##Distance2zero = []
        #for da in range(len(WtestP)):
            ##Distance2zero.append([pearsonr(WtestP[da], np.zeros(len(Wtest[da])))[0], euclidean(WtestP[da], np.zeros(len(Wtest[da])))])
            #for dq in range(len(WtrainP)):
                #SpredP[da,dq] = pearsonr(WtestP[da], WtrainP[dq])[0]
        #WtrainZ = (preprocessing.scale(WtrainZ.T)).T
        
        #SpredZ = np.zeros((len(WtrainZ),len(WtrainZ)))
        #for da in range(len(WtrainZ)):
            #for dq in range(len(WtrainZ)):
                #SpredZ[da,dq] = pearsonr(WtrainP[da], WtrainP[dq])[0]
                
        #SpredP2 = np.zeros((len(WtestP),len(WtrainP)))
        #for da in range(len(WtestP)):
            ##Distance2zero.append([pearsonr(WtestP[da], np.zeros(len(Wtest[da])))[0], euclidean(WtestP[da], np.zeros(len(Wtest[da])))])
            #for dq in range(len(WtrainP)):
                #SpredP2[da,dq] = pearsonr(SpredP[da], SpredZ[dq])[0]
        
        #SpredZ = np.zeros((len(WtestZ),len(WtrainZ)))
        #for da in range(len(WtestZ)):
            #for dq in range(len(WtrainZ)):
                #SpredZ[da,dq] = pearsonr(WtestZ[da], WtrainZ[dq])[0]
        #for sz, SZ in enumerate(SpredZ):
            #peaZ = pearsonr(SZ, SpredP[sz])[0]
            #print peaZ, pearsonr(SZ, SpredP2[sz])[0]
            ##if peaZ < 0.9:
                ##fig = plt.figure()
                ##ax = fig.add_subplot(111)
                ##ax.scatter(SpredP[sz], SpredZ[sz])
                ##plt.show()
        #sys.exit()


    else:
        print "Define similarity measure behind predS"
        sys.exit()
    
    if '--saveS' in sys.argv:
        print 'Similarity matrix saved'
        np.savetxt(os.path.splitext(Wname)[0]+'_'+smeasure+"-testset_S.txt", np.append(trainprots.reshape((1,len(trainprots))), Spred.astype(str), axis = 0).T, fmt="%s", header = ' '.join(testprots))
    
    if '--saveStest2testS' in sys.argv:
        if smeasure == 'dot': 
            Spredtest = np.dot(Ypred, Ytest.T)

        elif smeasure == 'pearson' or smeasure == 'pearsonreconstructed':
            Ypred = preprocessing.scale(Ypred)
            Spredtest = np.zeros((len(Wtest),len(Wtest)))
            for da in range(len(Ypred)):
                for dq in range(len(Ypred)):
                    Spredtest[da,dq] = pearsonr(Ypred[da], Ypred[dq])[0]
        
        elif smeasure == 'euclidean' or smeasure == 'euclideanreconstructed':
            Ypred = preprocessing.scale(Ypred)
            Spredtest = np.zeros((len(Ypred),len(Ytest)))
            for da in range(len(Ypred)):
                for dq in range(len(Ypred)):
                    Spredtest[da,dq] = euclidean(Ypred[da], Ypred[dq])

        elif smeasure == 'dotlatent': 
            Spredtest = np.dot(Wtest, Wtest.T)
        
        elif smeasure == 'pearsonlatent' or smeasure == 'pearsonlatentP':
            Wtest = (preprocessing.scale(Wtest.T)).T
            Spredtest = np.zeros((len(Wtest),len(Wtest)))
            for da in range(len(Wtest)):
                for dq in range(len(Wtest)):
                    Spredtest[da,dq] = pearsonr(Wtest[da], Wtest[dq])[0]
        
        elif smeasure == 'euclideanlatent' or smeasure == 'euclideanlatentP':
            Wtest = (preprocessing.scale(Wtest.T)).T
            Spredtest = np.zeros((len(Wtest),len(Wtest)))
            for da in range(len(Wtest)):
                for dq in range(len(Wtest)):
                    Spredtest[da,dq] = euclidean(Wtest[da], Wtest[dq])
            scaler = preprocessing.MinMaxScaler().fit(Spred.T)
            Spredtest = 1. - (scaler.transform(Spredtest.T)).T

        print 'Test similarity matrix saved'
        np.savetxt(os.path.splitext(Wname)[0]+'_'+smeasure+"-test2test_S.txt", np.append(testprots.reshape((1,len(testprots))), Spredtest.astype(str), axis = 0).T, fmt="%s", header = ' '.join(testprots))
    
    
    if '--predhighestcorr' in sys.argv:
        wtype = sys.argv[sys.argv.index("--predhighestcorr")+1]
        nfe = float(sys.argv[sys.argv.index("--predhighestcorr")+2])
        mixtype = sys.argv[sys.argv.index("--predhighestcorr")+3]
 
        if mixtype == 'top':
            def mixnum(s_test, ne):
                return int(ne)
        elif mixtype == 'cutoff':
            def mixnum(s_test, ne):
                neout = len(np.where(s_test >= ne)[0])
                return neout
 
        if wtype == "equal":
            def pred(s_test, y_train, ne):
                if ne == 0:
                    outvec = np.zeros(len(y_train[0]))
                else:
                    sort = np.argsort(-s_test)[:ne]
                    outvec = np.mean(y_train[sort], axis = 0)
                return outvec
        
        if wtype == "equalmax":
            def pred(s_test, y_train, ne):
                if ne == 0:
                    outvec = np.zeros(len(y_train[0]))
                else:
                    sort = np.argsort(-s_test)[:ne]
                    outvec = np.amax(y_train[sort], axis = 0)
                return outvec
            
        elif wtype == "weighted":
            def pred(s_test, y_train, ne):
                if ne == 0:
                    outvec = np.zeros(len(y_train[0]))
                else:
                    sort = np.argsort(-s_test)[:ne]
                    weight = s_test[sort]/np.sum(s_test[sort])
                    outvec = np.sum(y_train[sort].T*weight, axis = 1).T
                return outvec
        
        Ypred = []
        for i in range(len(testprots)):
            nge = mixnum(Spred[i], nfe)
            #print nge
            Ypred.append(pred(Spred[i], Ytrain, nge))
        Ypred = np.array(Ypred)
    
        if '--savetestcorrelation' in sys.argv:
            pearsont = np.zeros((len(testprots),2))
            for i in range(len(testprots)):
                pearsont[i,0] = pearsonr(Ypred[i], Ytest[i])[0]
                pearsont[i,1] = spearmanr(Ypred[i], Ytest[i])[0]
                print testprots[i], pearsont[i,0], pearsont[i,1]
            np.savetxt(os.path.splitext(Wname)[0]+'_'+smeasure+"-testset_reconstprofile"+wtype+str(nfe)+mixtype+"_pcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')
        
        if '--savetopintersection' in sys.argv:
            pearsont = np.zeros((len(testprots),2))
            perctop = int(len(Ypred[0])*0.01)
            for i in range(len(testprots)):
                pearsont[i,0] = len(np.intersect1d(np.argsort(Ypred[i])[-perctop:], np.argsort(Ytest[i])[-perctop:]))/float(perctop)
                pearsont[i,1] = len(np.intersect1d(np.argsort(Ypred[i])[-100:], np.argsort(Ytest[i])[-100:]))/(100.)       
                print testprots[i], pearsont[i,0], pearsont[i,1]
            np.savetxt(os.path.splitext(Wname)[0]+'_'+smeasure+"-testset_reconstprofile"+wtype+str(nfe)+mixtype+"_topintersection.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')
        
        if '--saveclosestcorrelation' in sys.argv:
                pearsont = []
                for i in range(len(testprots)):
                    pearsi = []
                    for j in range(len(trainprots)):
                        pearsi.append(pearsonr(Ypred[i], Ytrain[j])[0])
                    pearsont.append([trainprots[np.argmax(pearsi)], pearsi[np.argmax(pearsi)]] )
                #print pearsont
                np.savetxt(os.path.splitext(Wname)[0]+'_'+smeasure+wtype+str(nfe)+mixtype+"-reconstprofile_highcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), np.array(pearsont).astype(str), axis=1), fmt = '%s')
        
        if "--savetestprofiles" in sys.argv:
                np.savez_compressed(os.path.splitext(Wname)[0]+'_'+smeasure+"-testset_reconstprofile"+wtype+str(nfe)+mixtype+".npz", names=testprots, profiles=Ypred.T, kmers = Zfeatures)
        
        if "--savetestprofileslight" in sys.argv:
                tk = []
                for t, tprot in enumerate(testprots):
                    ksort = np.argsort(-Ypred[t])[:100]
                    tk.append(Zfeatures[ksort])
                np.savetxt(os.path.splitext(Wname)[0]+'_'+smeasure+"-testset_reconstprofile_light"+wtype+str(nfe)+mixtype+".txt", np.array(tk).T, fmt= '%s', header=' '.join(testprots))

        if '--savetopcorrelation' in sys.argv:
                    pearsont = np.zeros((len(testprots),2))
                    perctop = int(len(Ypred[0])*0.01)
                    for i in range(len(testprots)):
                        weight100 = np.zeros(len(Ypred[i]))
                        weight1p = np.zeros(len(Ypred[i]))
                        w1 = np.argsort(Ypred[i])[-perctop:]
                        w2 = np.argsort(Ytest[i])[-perctop:]
                        weight1p[w1] += 0.5
                        weight1p[w2] += 0.5
                        weight100[w1[-100:]] +=0.5
                        weight100[w2[-100:]] +=0.5
                        mask1p = weight1p > 0
                        mask100 = weight100 > 0
                        weight1p = weight1p[mask1p]
                        weight100 = weight100[mask100]
                        pearsont[i,0] = np.sum((Ypred[i][mask1p]-np.mean(Ypred[i][mask1p]))*weight1p*(Ytest[i][mask1p]-np.mean(Ytest[i][mask1p])))/(np.sum(weight1p)*np.std(Ypred[i][mask1p])*np.std(Ytest[i][mask1p]))
                        pearsont[i,1] = np.sum((Ypred[i][mask100]-np.mean(Ypred[i][mask100]))*weight100*(Ytest[i][mask100]-np.mean(Ytest[i][mask100])))/(np.sum(weight100)*np.std(Ypred[i][mask100])*np.std(Ytest[i][mask100]))
                    np.savetxt(os.path.splitext(Wname)[0]+'_'+smeasure+"-testset_reconstprofile"+wtype+str(nfe)+mixtype+"_topcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')
                
                
                
                
                
                
                
