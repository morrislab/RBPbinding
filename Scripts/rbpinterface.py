#rpbinterface.py
import numpy as np
import sys, os
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from itertools import imap
from operator import eq
from joblib import Parallel, delayed
import multiprocessing
from scipy.sparse import csc_matrix, issparse, hstack
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62



aa = ['A', 'R', 'N', 'D', 'P', 'V', 'I', 'C', 'Y', 'H', 'T', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S']


def readinmask(maskfile):
    obj = open(maskfile, 'r').readlines()
    pseq = []
    pname = []
    pinter = []
    pinterdist = []
    interrna = []
    pinterback = []
    pinterresidue = []
    for l, line in enumerate(obj):
        if line[0] == '>':
            pname.append(line.strip()[1:])
            pcseq = np.array(obj[l+1].strip().split(','))
            pcinter = np.array(obj[l+2].strip().split(','), dtype = int)
            cinterrna = obj[l+5].strip().split(',')
            pcinterback = obj[l+4].strip().split(',')
            pcinterresidue = obj[l+3].strip().split(',')
            pcinterdist = np.array(obj[l+6].strip().split(','), dtype = float)
            
            pcmask = np.concatenate([[-1],np.where(~np.isin(pcseq, aa))[0], [len(pcseq)]])
            qseq = []
            qinter = []
            qrna = []
            qback = []
            qres = []
            qdist = []
            for q in range(len(pcmask)-1):
                qseq.append(''.join(pcseq[pcmask[q]+1:pcmask[q+1]]))
                qinter.append(pcinter[pcmask[q]+1:pcmask[q+1]])
                qrna.append(cinterrna[pcmask[q]+1:pcmask[q+1]])
                qback.append(pcinterback[pcmask[q]+1:pcmask[q+1]])
                qres.append(pcinterresidue[pcmask[q]+1:pcmask[q+1]])
                qdist.append(pcinterdist[pcmask[q]+1:pcmask[q+1]])
            pseq.append(qseq)
            pinter.append(qinter)
            interrna.append(qrna)
            pinterback.append(qback)
            pinterresidue.append(qres)
            pinterdist.append(qdist)
    return pseq, pname, pinter, pinterdist, interrna, pinterback, pinterresidue
            
     

def readinfasta(fastafile):
    obj = open(fastafile, 'r').readlines()
    pseq = []
    pname = []
    for l, line in enumerate():
        if line[0] == '>':
            pname.append(line.strip()[1:])
            pseq.append(''.join(np.array(obj[l+1].strip().split(','))))
    return pseq, pname



   

def add_sing(arr, sing):
    outarr = []
    for ar in arr:
        for si in sing:
            outarr.append(ar+si)
    return outarr

def createkmer(sequence, klen, gapsize, ktype, allk = False):
    if int(allk)> 1:
        start = int(allk)
    elif allk == True:
        start = 1
    else:
        start = klen
    print ktype, klen, gapsize
    if ktype == 'gapped':
        cent = (klen - gapsize)/2
        centstart = []
        for cs in range(cent, klen -cent):
            centstart.append(cs)
            if cs +gapsize == klen -cent:
                break
        kmers = ['']
        for i in range(klen - gapsize):
            kmers = add_sing(kmers, aa)
        kmers = np.array(kmers)
        #print kmers
        featmat = []
        for i, seqarr in enumerate(sequence):
            #print i
            farray = []
            for seq in seqarr:
                row = [] 
                col = [] 
                data = []
                shape = (len(kmers), len(seq))
                for j in range(len(seq)+1-klen):
                    fkmer = seq[j:j+klen]
                    for cent in centstart:
                        row.append(list(kmers).index(fkmer[:cent]+fkmer[cent+gapsize:]))
                        col.append(j)
                        data.append(1)
                featarray = csc_matrix((data, (row, col)), shape=shape)
                farray.append(featarray.toarray())
            #if len(seqarr) > 1:
                #farray = np.concatenate(farray, axis = -1)
            #else:
                #farray = farray[0]
            featmat.append(farray)
        okmers = []
        for k, kmer in enumerate(kmers):
            okmers.append('G'+str(gapsize)+kmer)
        kmers = okmers
    elif ktype == 'mismatch':
        kmers = ['']
        for i in range(klen - gapsize):
            kmers = add_sing(kmers, aa)
        kmers = np.array(kmers)
        
        # define the position which are kept
        def mix(inseq, kma):
            och = True
            if inseq[-1] != kma:
                inseq[-1]+=1
            else:
                if len(inseq) > 1:
                    out, och = mix(inseq[:-1], kma-1)
                    inseq[-1] = inseq[-2] +1
                else:
                    och = False
            return inseq, och
        
        listpositions = []
        allpos = np.arange(klen)
        posar = np.arange(gapsize)
        ocheck = True
        while ocheck:
            listpositions.append(np.setdiff1d(allpos,posar))
            posar, ocheck = mix(posar, klen-2) # Never get cover the last kmer because otherwise these kmers are counted twice
        listpositions = np.array(listpositions)    
        
        featmat = [] 
        for i, seqarr in enumerate(sequence):
            farray = []
            for seq in seqarr:
                row = [] 
                col = [] 
                data = []
                shape = (len(kmers), len(seq))
                for j in range(len(seq)+1-klen):
                    fkmer = list(seq[j:j+klen])
                    for lpos in listpositions:
                        mikmer = ''.join(np.array(list(fkmer))[lpos])
                        row.append(list(kmers).index(mikmer))
                        col.append(j)
                        data.append(1)
                featarray = csc_matrix((data, (row, col)), shape=shape)
                farray.append(featarray.toarray())
            featmat.append(farray)
        okmers = []
        for k, kmer in enumerate(kmers):
            okmers.append('M'+str(gapsize)+kmer)
        kmers = okmers            
    elif ktype == 'posmismatch':
        kmerfeat = []
        for kend in range(start-gapsize,k+1-gapsize):    
            kmers = ['']
            for i in range(kend):
                kmers = add_sing(kmers, aa)
            kmerfeat.append(kmers)
        print 'Generate kmer'
        kmers = np.concatenate(kmerfeat)
        rkmerfeat = np.copy(kmers)
        kmers = []
        for g in range(gapsize):
            gappedkmerfeat = []
            for kmer in rkmerfeat:
                gkmer = []
                for i in range(1,len(kmer)):
                    gkmer.append(kmer[:i]+'X'+kmer[i:])
                gappedkmerfeat.append(gkmer)
            rkmerfeat = np.unique(np.concatenate(gappedkmerfeat))
        kmers = rkmerfeat
        print 'kmers done ...'

        def findkmer(seqar, kmerset, pomask, iy):
            koutarray = []
            ml = -1
            for l in range(start,klen+1):
                ml += 1
                masks = pomask[ml]
                for seq in seqar:
                    row = [] 
                    col = [] 
                    data = []
                    shape = (len(kmers), len(seq))
                    for j in range(len(seq)+1-l):
                        ckm = np.array(list(seq[p:p+l]))
                        #print ckm
                        for mask in masks:
                            cckm = np.copy(ckm)
                            cckm[mask] = 'X'
                            #print cckm
                            row.append(list(kmerset).index(''.join(cckm)))
                            col.append(j)
                            data.append(1)
                    featarray = csc_matrix((data, (row, col)), shape=shape)
                    koutarray.append(featarray.toarray())
            return koutarray, iy

        def makemask(start, k, gapsize):
            omasks = []
            for l in range(start,k+1):
                masks = []
                for h in range(1,l-1):
                    mask = np.zeros(l) > 0.
                    mask[h] = True
                    masks.append(mask)
                rmasks = np.copy(masks)
                for g in range(gapsize-1):
                    gmasks = []
                    for mask in rmasks:
                        for h in range(1,l-1):
                            cmask = np.copy(mask)
                            cmask[h] = True
                            if np.sum(cmask) == g + 2:
                                gmasks.append(cmask)
                    rmasks = np.unique(gmasks, axis = 0)
                omasks.append(rmasks)
            return omasks
        
        featmat = []
        tmasks = makemask(start,klen,gapsize)
        print 'Get k-mers...'
        if mprocessing:
            results, order = Parallel(n_jobs=num_cores)(delayed(findkmer)(sequence[i], kmers, tmasks, i) for i in range(len(sequence)))
            print 'Distribute multiprocessing results'
            order = np.argsort(order)
            for rpos in order:
                featmat.append(results[rpos])
        else:
            for i, seqarr in enumerate(sequence):
                res, r = findkmer(seqarr, kmers, tmasks, i)
                featmat.append(res)
    else:    
        kmerfeat = []
        for kend in range(start,klen+1):    
            kmers = ['']
            for i in range(kend):
                kmers = add_sing(kmers, aa)
            kmerfeat.append(kmers)
        kmers = np.concatenate(kmerfeat)
        #### search for kmers
        featmat = []
        for i, seqarr in enumerate(sequence):
            #print i
            farray = []
            for s, seq in enumerate(seqarr):
                row = [] 
                col = [] 
                data = []
                shape = (len(kmers), len(seq))
                for kend in range(start,klen+1): 
                    for j in range(len(seq)+1-kend):
                        row.append(list(kmers).index(seq[j:j+kend]))
                        col.append(j)
                        data.append(1)
                featarray = csc_matrix((data, (row, col)), shape=shape)
                #print featarray
                farray.append(featarray.toarray())
            #if len(seqarr) > 1:
                #farray = np.concatenate(farray, axis = -1)
            #else:
                #farray = farray[0]
            featmat.append(farray)
    #print featmat
    return kmers, featmat


# weighting can be divided by the number of hits per sequence, no weighting
class classifier():
    def __init__(self, classtype, parameters):
        # introduce weighting
        self.classtype = classtype
        self.parameters = parameters
        if classtype == 'RF':
            if len(parameters) != 0:
                self.nesti = int(parameters[0])
            else:
                self.nesti = 100
            self.cl = RandomForestClassifier(n_estimators=self.nesti, n_jobs=num_cores)
        elif classtype == 'logit':
            if len(parameters) != 0:
                self.penalty = parameters[0]
                self.tol = float(parameters[1])
                self.cc = float(parameters[2])
                self.classweight = parameters[3] # balanced
                self.solver = parameters[4]
                self.maxiter = int(parameter[5])
            else:
                self.penalty = 'l2'
                self.tol = 0.0000001
                self.cc = 1000000. # inverse of regularization strength
                self.classweight = None # balanced
                self.solver = 'saga'
                self.maxiter = 10000
            self.cl = LogisticRegression(penalty=self.penalty, tol=self.tol, C=self.cc, fit_intercept=True, solver=self.solver, max_iter=self.maxiter, n_jobs=num_cores)
        elif classtype == 'svm':
            if len(parameters) != 0:
                self.kernel = parameters[0]
                self.cc = float(parameters[1])
                self.tol = float(parameters[2])
                self.cache = int(parameters[3])
            else:
                self.kernel = 'linear'
                self.cc = 1000.0
                self.tol = 0.001
                self.cache = 8000
            self.cl = SVC(C=self.cc, kernel=self.kernel, probability = True, tol=self.tol, cache_size=self.cache, max_iter=-1, gamma='auto')
        elif classtype == 'GB':
            if len(parameters) != 0:
                self.loss = parameters[0]
                self.lrate = float(parameters[1])
                self.nesti = int(parameters[2])
                self.maxdepth = int(parameters[3])
                self.tol = float(parameters[4])
            else:
                self.loss = 'deviance'
                self.lrate = 0.1
                self.nesti = 100
                self.maxdepth = 3
                self.tol = 0.0001
            self.cl = GradientBoostingClassifier(loss=self.loss, learning_rate=self.lrate, n_estimators=self.nesti, subsample=1.0, min_samples_split=2, min_samples_leaf=1, max_depth=self.maxdepth, tol=self.tol)
        print 'Fit data', classtype, parameters
    def train(self, xtr, ytr, dweights):
        #print xtr, ytr
        if dweights is None:
            self.cl.fit(xtr,ytr)
        else:
            self.cl.fit(xtr, ytr, dweights)
    def predict(self, xte):
        outclass = self.cl.predict(xte)
        outscore = self.cl.predict_proba(xte)
        return outclass, outscore[:,-1]
        
    def coef(self):
        if self.classtype == 'RF':
            fweight = self.cl.feature_importances_
        if self.classtype == 'logit':
            fweight = self.cl.coef_[0]
        if self.classtype == 'svm':
            if self.kernel == 'linear':
                ### Need to optimize for rbf kernel, only applicable to linear kernel
                fweight = self.cl.coef_[0]
            else:
                fweight = np.zeros(len(self.cl.support_vectors_[0]))
        if self.classtype == 'GB':
            fweight = self.cl.feature_importances_
        return fweight
        
    
    
    

def makedataset(seqrep, yin, window, combine, normalizex = False):
    xdata = []
    ydata = []
    
    for i, yi in enumerate(yin):
        ydata.append(np.concatenate(yi))
    
    if combine == 'sum':
        def combf(inarr, spars):
            if spars:
                return csc_matrix(inarr.sum(axis = 1))
            else:
                return np.sum(inarr, axis = 1)
            
    elif combine == 'concatenate':
        def combf(inarr, spars):
            if spars:
                spashape = inarr.shape
                return inarr.reshape(spashape[0] * spashape[1])
            else:
                return inarr.flatten()
    
    for s, seq in enumerate(seqrep):
        nseq = []
        #print s
        if issparse(seq[0]):
            for z, sez in enumerate(seq):
                #print z
                sezshape = sez.shape
                sezin = hstack([csc_matrix((sezshape[0], int(window/2.))), sez, csc_matrix((sezshape[0], window))])
                for p in range(sezshape[-1]):
                    xseq = combf(sezin[:,p:p+window], True)
                    if normalizex:
                        xseq = xseq/xseq.sum()
                    #print xseq.shape
                    nseq.append(xseq)
        else:
            for z, sez in enumerate(seq):
                sezshape = np.shape(sez)
                sezin = np.concatenate([np.zeros((sezshape[0], int(window/2.))), sez, np.zeros((sezshape[0], window))], axis = 1)
                for p in range(sezshape[-1]):
                    xseq = combf(sezin[:,p:p+window], False)
                    if normalizex:
                        xseq = xseq/np.sum(xseq)
                        xseq[np.isinf(xseq)] = 0.
                        xseq[np.isnan(xseq)] = 0.
                    nseq.append(xseq)
        #print np.shape(nseq), len(ydata[s])
        xdata.append(nseq)
    #sys.exit()
    return xdata, ydata

def seqweight(ytr, protweight):
    weight = []
    for t, yt in enumerate(ytr):
        wyt = np.ones(len(yt))
        if protweight == 'protlen':
            wyt = wyt/float(len(wyt))
        if protweight == 'numres':
            wyt = wyt/(5.+np.sum(wyt))    
        if protweight == 'protlen-numres':
            wyt = wyt/float(len(wyt))
            wyt = wyt/(5.+np.sum(wyt))
        weight.append(wyt)
    return weight

    
def interfacefile(oname, ypred, tseqs, tprots):
    wobj = open(oname+'interface.fasta', 'w')
    for t, tprot in enumerate(tprots):
        wobj.write('>'+tprot+'\n'+''.join(np.array(tseqs[t]))+'\n'+','.join(ypred[t][0].astype(str))+'\n'+','.join(np.around(ypred[t][1],3).astype(str))+'\n')
    wobj.close()
    
def precisionrecall(y_score, y_test):
    from sklearn.metrics import average_precision_score
    praverage = []
    for yi in range(len(y_test)):
        average_precision = average_precision_score(y_test[yi], y_score[yi][1])
        praverage.append(average_precision)
    return praverage

def roc(y_score, y_test):
    from sklearn.metrics import roc_curve, auc
    aucaverage = []
    for yi in range(len(y_test)):
        fpr, tpr, _ = roc_curve(y_test[yi], y_score[yi][1])
        roc_auc = auc(fpr, tpr)
        aucaverage.append(roc_auc)
    return aucaverage
    
def saveperformance(oname, tprots, pr, auc):
    wobj = open(oname+'-prauc.dat', 'w')
    wobj.write('# RBP PR AUC\n')
    for o, on in enumerate(tprots):
        wobj.write(on+' '+ str(pr[o])+' '+str(auc[o])+'\n')
    wobj.close()
    
    
if __name__ == '__main__':
    
    if '--mprocessing' in sys.argv:
        mprocessing = True
        num_cores = int(sys.argv[sys.argv.index('--mprocessing')+1])
    else:
        mprocessing = False
        num_cores = 1
    
    
    
    if '--trainmodel' in sys.argv:
        protseq, protname, protinter, protinterdist, interfacerna, protinterback, protinterresidue = readinmask(sys.argv[sys.argv.index('--trainmodel')+1])
        train = True
        #print protseq, protname, protinter, protinterdist, interfacerna, protinterback, protinterresidue
    elif '--makepredictions' in sys.argv:
        protseq, protname = readinfasta(sys.argv[sys.argv.index('--makepredictions')+1])
        train = False
    
    if '--outname' in sys.argv:
        outself = True
        outname = sys.argv[sys.argv.index('--outname')+1]
    else:
        outself = False
        outname = 'Interfacepred'
    
    ## instead of k-mers, could also use different physico-chemical properties
    ## or could use different clusters and for each cluster have one-hot encoding. Feature vector is concatenated featurevector of all surrounding residues
    if '--kmer' in sys.argv:
        kmerlen = int(sys.argv[sys.argv.index('--kmer')+1])
        gaptype = sys.argv[sys.argv.index('--kmer')+2] ### all k<=kmerlen k-mers can be determined, allow blosum smoothing
        gaplen = int(sys.argv[sys.argv.index('--kmer')+3])
        outname += '-k'+str(kmerlen)+gaptype+str(gaplen)
        if '--allk' in sys.argv:
            allk = True
            outname += 'all'
        else:
            allk = False
        kmers, kmermat = createkmer(protseq, kmerlen, gaplen, gaptype, allk)

    elif '--kmerfile' in sys.argv:
        # should also be the option to provide precaluculated file
        # precaluculated files can contain k-mers from orthologs
        # create new scripts that do read out k-mer per position
        # create script that generates feature-file from pssm, singleAA, asapred, disorderpred, secstrucpred
        kfile = np.load(sys.argv[sys.argv.index('--kmerfile')+1],allow_pickle=True)
        print kfile.files
        kmers = kfile['features']
        kmermatin = kfile['featmat']
        kprotnames = kfile['protnames']
        ksequences = kfile['sequences']
        # rearrange to input data and check if kmer-profiles possess correct length
        sort = []
        if outself == False:
            outname += '-'+os.path.splitext(os.path.split(sys.argv[sys.argv.index('--kmerfile')+1])[1])[0]
        kmermat = []
        # sometimes the masks don't fit the whole domain because X have been removed, and so they need to be aligned 
        for p, prot in enumerate(protname):
            kind = list(kprotnames).index(prot)
            cprotseq = ''.join(protseq[p])
            kfeatmat = kmermatin[kind].T
            if ksequences[kind] != cprotseq:
                print len(ksequences[kind]), len(cprotseq)
                print np.shape(kfeatmat)
                alignment = pairwise2.align.globalds(ksequences[kind], cprotseq, matrix, -11, -1)
                alignmask = np.array(list(alignment[0][1])) != '-'
                kfeatmat = kfeatmat[:,alignmask]
            sort.append(kind)
            kmermat.append([kfeatmat])
        kprotnames = kprotnames[sort]
    
    else:
        kmerlen = 2
        gaptype = None
        gaplen = 0
        kmers, kmermat = createkmer(protseq, kmerlen, gaptype, gaplen, False)
        outname += 'k'+str(kmerlen)+gaptype+str(gaplen)
    
    # second k-mertype can be added
    if '--addkmer' in sys.argv:
        kmerlen = int(sys.argv[sys.argv.index('--addkmer')+1])
        gaptype = sys.argv[sys.argv.index('--addkmer')+2] ### all k<=kmerlen k-mers can be determined, allow blosum smoothing
        gaplen = int(sys.argv[sys.argv.index('--addkmer')+3])
        if '--addallk' in sys.argv:
            allk = True
        else:
            allk = False
        akmers, akmermat = createkmer(protseq, kmerlen, gaplen, gaptype, allk)
        #print akmers, akmermat
        outname += '-addk'+str(kmerlen)+gaptype+str(gaplen)
    elif '--addkmerfile' in sys.argv:
        kfile = np.load(sys.argv[sys.argv.index('--addkmerfile')+1],allow_pickle=True)
        print kfile.files
        akmers = kfile['features']
        akmermatin = kfile['featmat']
        akprotnames = kfile['protnames']
        aksequences = kfile['sequences']
        # rearrange to input data and check if kmer-profiles possess correct length
        sort = []
        if outself == False:
            outname += '-'+os.path.split(sys.argv[sys.argv.index('--addkmerfile')+1])[1].split('.')[0]
        akmermat = []
        for p, prot in enumerate(protname):
            kind = list(akprotnames).index(prot)
            cprotseq = ''.join(protseq[p])
            akfeatmat = akmermatin[kind].T
            if aksequences[kind] != cprotseq:
                print len(aksequences[kind]), len(cprotseq)
                alignment = pairwise2.align.globalds(aksequences[kind], cprotseq, matrix, -11, -1)
                alignmask = np.array(list(alignment[0][1])) != '-'
                akfeatmat = akfeatmat[:,alignmask]
            sort.append(kind)
            akmermat.append(np.array([akfeatmat]))
        akprotnames = akprotnames[sort]
        
    if '--addkmer' in sys.argv or '--addkmerfile' in sys.argv:
        kmers = np.append(kmers, akmers)
        #print kmers
        for p in range(len(protname)):
            for s in range(len(kmermat[p])):
                kmermat[p][s] = np.append(kmermat[p][s], akmermat[p][s], axis = 0)
            
    if '--trainingset' in sys.argv:
        # either give file that determines which sequences are used for training or use percentage of permutation
        numset = sys.argv[sys.argv.index("--trainingset")+1]
        trainfile = sys.argv[sys.argv.index("--trainingset")+2]
        if outself == False:
            outname += '-TRain'+os.path.splitext(os.path.split(trainfile)[1])[0]+numset
        else:
            outname += 'TRain'+numset
        print "trainset:", trainfile, numset
        tobj = open(trainfile, 'r')
        tlines = tobj.readlines()
        for i, tline in enumerate(tlines):
            if tline[:6]=="###Set" and tline.strip().split()[-1] == numset:
                if tlines[i+1][:7] == '##Train':
                    trainprots = np.array(tlines[i+2].strip().split())
                if tlines[i+3][:6] == '##Test':
                    testprots = np.array(tlines[i+4].strip().split())
                if tlines[i+4][:6] == '##Test':
                    testprots = np.array(tlines[i+5].strip().split())
                break

        if np.sum(np.isin(protname, trainprots)) == 0:
            for t, trainp in enumerate(trainprots):
                trainprots[t] = trainp[1:]
            for t, trainp in enumerate(testprots):
                testprots[t] = trainp[1:]
        if np.sum(np.isin(protname, trainprots)) == 0:
            print 'trainprots not in protnames', trainprots, testprots, protname
            sys.exit()
        ### sort into training and test set
        indexestrain = []
        indexestest = []
        for i, ina in enumerate(protname):
            if ina in trainprots:
                indexestrain.append(i)
            if ina in testprots:
                indexestest.append(i)
        indexestrain = np.array(indexestrain)
        indexestest = np.array(indexestest)    
        
        trainprots = np.array(protname)[indexestrain]
        testprots = np.array(protname)[indexestest]
        
        trainseqs = [] 
        testseqs = []
        traininter = []
        testinter = []
        feattrain = [] 
        feattest = []
        for i in indexestrain:
            feattrain.append(kmermat[i])
            trainseqs.append(protseq[i])
            traininter.append(protinter[i])
        for i in indexestest:
            feattest.append(kmermat[i])
            testseqs.append(protseq[i])
            testinter.append(protinter[i])
            
    elif '--crossvalidation' in sys.argv:
        cross = int(sys.argv[sys.argv.index('--crossvalidation')+1])
        crosset = int(sys.argv[sys.argv.index('--crossvalidation')+2])
        outname += '-cv'+str(crosset)+'-'+str(cross)
        lentest = len(protname)/cross
        indexestrain = np.arange(len(protname))
        indexestest = np.arange(len(protname)-crosset*lentest,len(protname)-crosset*lentest+lentest)
        indexestrain = np.delete(indexestrain, indexestest)
        #print indexestrain, indexestest
        trainprots = np.array(protname)[indexestrain]
        testprots = np.array(protname)[indexestest]
        
        trainseqs = [] 
        testseqs = []
        traininter = []
        testinter = []
        feattrain = [] 
        feattest = []
        for i in indexestrain:
            feattrain.append(kmermat[i])
            trainseqs.append(protseq[i])
            traininter.append(protinter[i])
        for i in indexestest:
            feattest.append(kmermat[i])
            testseqs.append(protseq[i])
            testinter.append(protinter[i])
        
    ### define model: svm, randomforest, log regression, gradient boost
    if '--classifier' in sys.argv:
        model = sys.argv[sys.argv.index('--classifier')+1] # Classifier
        windowsize = int(sys.argv[sys.argv.index('--classifier')+2]) # sequence window needs to uneven
        featuremix = sys.argv[sys.argv.index('--classifier')+3] # sum or concatenate features in window
        protweight = sys.argv[sys.argv.index('--classifier')+4] # weight sequence by length or number of positive
        xnorm = bool(sys.argv[sys.argv.index('--classifier')+5]) # useless because can be done before giving to model
        if xnorm:
            xnoou = 'xnorm'
        else:
            xnoou = ''
        parameters = [] 
        outname += '-'+model+'-wz'+str(windowsize)+featuremix+'-'+protweight+xnoou
        print 'Generate trainfeatures', windowsize, featuremix, xnorm
        xtrain, ytrain = makedataset(feattrain, traininter, windowsize, featuremix, xnorm)
        xtest, ytest = makedataset(feattest, testinter, windowsize, featuremix, xnorm)
    if train:
        if protweight == 'protlen-numres' or protweight == 'protlen' or protweight == 'numres':
            pdweight = seqweight(ytrain, protweight)
            pdweight = np.concatenate(pdweight, axis = 0)
        else:
            pdweight = None
        ### perform training
        if '--modelparams' in sys.argv:
            parameters = sys.argv[sys.argv.index('--modelparams')+1].split(',')
        cmodel = classifier(model, parameters)
        cmodel.train(np.concatenate(xtrain, axis = 0), np.concatenate(ytrain, axis = 0), pdweight)
        featureweights = cmodel.coef()
        #print len(featureweights)
        fsort = np.argsort(-featureweights)
        if featuremix == 'sum':
            featureout = kmers
        elif featuremix == 'concatenate':
            featureout = []
            for i in range(windowsize):
                for kmer in kmers:
                    featureout.append(str(i)+"#"+kmer)
        
        np.savetxt(outname +'-fweight.dat', np.append([np.array(featureout)[fsort]],[featureweights[fsort]], axis = 0).T.astype(str), fmt = '%s')
        
    else:
        loadclassifier = True
    
    #ypredtrain = []
    #for xt in xtrain:
        #ypredtrain.append(cmodel.predict(xt))
    ypredtest = []
    for j, xt in enumerate(xtest):
        ypredtest.append(cmodel.predict(xt))
        print j, len(xtest)
    
    # create fasta file for prediction scores and binary form
    interfacefile(outname+'-test', ypredtest, testseqs, testprots)
    #interfacefile(outname+'-intfacetrain', ypredtrain, trainseqs, trainprots)
    # determine auc, prauc for each test protein and save both
    
    prtest = precisionrecall(ypredtest, ytest)
    auctest = roc(ypredtest, ytest)
    
    saveperformance(outname, testprots, prtest, auctest)
    
    print 'AUC', np.mean(auctest), '+-', np.std(auctest)
    print 'PRAUC', np.mean(prtest), '+-', np.std(prtest)
        
        
        
        
