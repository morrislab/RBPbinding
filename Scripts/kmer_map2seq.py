#mapkmerscore2sequence.py
# map significant kmers to the protein sequence
import numpy as np
import scipy
import sys, os
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62
from Bio import pairwise2
import operator
from joblib import Parallel, delayed
import multiprocessing
import matplotlib.pyplot as plt

def map2seq(kweights, pkmer, pseq, ly):
    print ly
    posweights = np.zeros(len(pseq))
    print pseq
    if '*' in pseq:
        seqs = pseq.split('*')
    else:
        seqs = [pseq]
    print len(pseq), seqs, len(posweights)
    for pp, pep in enumerate(pkmer):
        if kweights[pp] != 0:
            if 'X' in pep:
                kmask = np.array(list(pep)) != 'X'
                klen = len(pep)
                kmcomp = np.array(list(pep))[kmask]
                fk = 0.
                reweight = []
                offset = 0
                for q, qseq in enumerate(seqs): 
                    for pr in range(len(qseq) - klen +1):
                        compseq = np.array(list(qseq[pr:pr+klen]))
                        compseq = compseq[kmask]
                        
                        if np.array_equal(compseq, kmcomp):
                            fk += 1.
                            weloc = np.arange(klen) + offset + pr
                            reweight.append(weloc[kmask])
                    offset += len(qseq) + 1
            else:
                klen = len(pep)
                fk = 0.
                reweight = []
                offset = 0
                for q, qseq in enumerate(seqs):
                    for pr in range(len(qseq) - klen +1):
                        compseq = qseq[pr:pr+klen]
                        if compseq == pep:
                            fk += 1.
                            weloc = np.arange(klen) + offset + pr
                            reweight.append(weloc)
                    offset += len(qseq) + 1
                
            if len(reweight) > 0:
                #print reweight, pep, pseq
                reweight = np.concatenate(reweight)
                #print reweight
                #print pep
                #print np.array(list(pseq))[reweight] 
                posweights[reweight] += kweights[pp]/fk

    return ly, posweights



if '--kmerfeatures' in sys.argv:
    common_name = sys.argv[sys.argv.index('--kmerfeatures')+1]
    #Ppred = np.load(common_name)['Ppred']
    #pnames = np.load(common_name)['pnames']
    #Pfeatures =np.load(common_name)['Pfeatures']
    Ppred = np.load(common_name)['profiles']
    shp = np.shape(Ppred)
    if shp[0] > shp[1]:
        Ppred = Ppred.T
    print np.shape(Ppred)
    pnames = np.array(np.load(common_name)['names'])
    Pfeatures =np.load(common_name)['kmers']
    
if '--minlen' in sys.argv:
    minlen = sys.argv[sys.argv.index('--minlen')+1]
    common_name = os.path.splitext(common_name)[0]+'-mlfeat'+minlen+os.path.splitext(common_name)[1]
    mylen = np.vectorize(len)
    Pfeaturelen = mylen(Pfeatures)
    keep = Pfeaturelen >= int(minlen)
    Pfeatures = Pfeatures[keep]
    Ppred = Ppred[:,keep]
    
if '--normalize' in sys.argv:
    common_name = os.path.splitext(common_name)[0]+'-norm'+os.path.splitext(common_name)[1]
    print 'normalize'
    mylen = np.vectorize(len)
    Pfeaturelen = mylen(Pfeatures)
    uniquelen = np.unique(Pfeaturelen)
    if len(uniquelen) > 1:
        for ul, unilen in enumerate(uniquelen):
            pmask = Pfeaturelen == unilen
            mean = np.mean(Ppred[:, pmask])
            std = np.std(Ppred[:, pmask])
            Ppred[:, pmask] -= mean
            Ppred[:, pmask] /= std
    else:
        mean = np.mean(Ppred)
        std = np.std(Ppred)
        print mean, std
        #print Ppred[0,0]
        Ppred -= np.mean(Ppred)
        Ppred /= np.std(Ppred)
        print np.amin(Ppred), np.amax(Ppred)
        print np.amin(np.amax(Ppred, axis = 1)), np.amax(np.amin(Ppred, axis = 1))
    #sys.exit()

kweight = np.copy(Ppred)


if '--original_kmerfeatures' in sys.argv:
    orig_name = sys.argv[sys.argv.index('--original_kmerfeatures')+1]
    orobj = np.load(orig_name) 
    Porig = orobj['features']
    ornames = np.array(orobj['expnames'])
    orkmers = orobj['kmers']
    # Only keep the kmer weights of k-mers that occur in the sequence
    # alignfeatures to each other
    psortrec = np.argsort(Pfeatures)
    psortorig = np.argsort(orkmers)
    Ppred = Ppred[:, psortrec]
    Porig = Porig[psortorig]
    Pfeatures = Pfeatures[psortrec]
    orkmers = orkmers[psortorig]
    preckeep = np.isin(Pfeatures, orkmers)
    porkeep = np.isin(orkmers, Pfeatures)
    Ppred = Ppred[:, preckeep]
    Porig = Porig[porkeep]
    Pfeatures = Pfeatures[preckeep]
    orkmers = orkmers[porkeep]
    print 'Filtering kmers for original set', np.array_equal(Pfeatures, orkmers)
    
    print 'Found Proteins', len(pnames), len(np.intersect1d(ornames, pnames))
    mask = np.isin(pnames, ornames)
    pnames = pnames[mask]
    Ppred = Ppred[mask]
    kweight = kweight[mask]
    for p, pname in enumerate(pnames):
        op = list(ornames).index(pname)
        Ppred[p] = Ppred[p] * (Porig[:, op] > 0).astype(float)
        


if '--domaininfo' in sys.argv:
    infofile = sys.argv[sys.argv.index('--domaininfo')+1]
    startout = os.path.splitext(os.path.split(infofile)[1])[0]
    sobj = open(infofile, 'r')
    protnames = []
    expnames = []
    sequences = []
    for i, sline in enumerate(sobj):
        sline = sline.strip().split()
        protnames.append(sline[0].split("||")[0])
        expnames.append(sline[0].split("||")[1])
        if '_' in protnames[0]:
            for pm, protname in enumerate(protnames):
                protnames[pm] =protname.split('_')[0]
        ndom = int(sline[1])
        sequences.append('*'.join(sline[-ndom:]))
    keep = []
    pnamekeep = []
    for p, pname in enumerate(pnames):
        if pname in expnames:
            keep.append(expnames.index(pname))
            pnamekeep.append(p)
    sequences = np.array(sequences)[keep]
    kweight = kweight[pnamekeep]
    pnames = pnames[pnamekeep]
            
elif '--sequencefasta' in sys.argv:
    sequencefile = sys.argv[sys.argv.index("--sequencefasta")+1]
    startout = os.path.splitext(os.path.split(sequencefile)[1])[0]
    sobj = open(sequencefile, 'r').readlines()
    protnames = []
    expnames = []
    sequences = []
    for i, sline in enumerate(sobj):
            sline = sline.strip()
            if sline[0] == '>':
                if '||' in sline:
                    protnames.append(sline.split("||")[0][1:])
                    expnames.append(sline.split("||")[1])
                else:
                    protnames.append(sline[1:])
                    expnames.append(sline[1:])
                objseq = sobj[i+1].strip()
                #if '*' in objseq:
                    #objseq = objseq.split('*')
                #else:
                    #objseq = [objseq]
                sequences.append(objseq)
    keep = []
    pnamekeep = []
    for p, pname in enumerate(pnames):
        if pname in expnames:
            keep.append(expnames.index(pname))
            pnamekeep.append(p)
    sequences = np.array(sequences)[keep]
    kweight = kweight[pnamekeep]
    pnames = pnames[pnamekeep]
   
elif '--domainfasta' in sys.argv:
    sequencefile = sys.argv[sys.argv.index("--domainfasta")+1]
    startout = os.path.splitext(os.path.split(sequencefile)[1])[0]
    sobj = open(sequencefile, 'r').readlines()
    dprotnames = []
    dexpnames = []
    domainnames = []
    dsequences = []
    for i, sline in enumerate(sobj):
        sline = sline.strip()
        if sline[0] == '>':
            domainnames.append(sline.split('__')[-1])
            sline = sline.split('__')[0]
            if '||' in sline:
                dprotnames.append(sline.split("||")[0][1:])
                dexpnames.append(sline.split("||")[1])
            else:
                dprotnames.append(sline[1:])
                dexpnames.append(sline[1:])
            objseq = sobj[i+1].strip()
            dsequences.append(objseq)
    dexpnames = np.array(dexpnames)
    dsequences = np.array(dsequences)
    expnames = np.unique(dexpnames)
    sequences = []
    protnames = []
    for e, expn in enumerate(expnames):
        sind = np.where(dexpnames == expn)[0]
        sequences.append('*'.join(dsequences[sind]))
        protnames.append(dprotnames[sind[0]])
    keep = []
    pnamekeep = []
    for p, pname in enumerate(pnames):
        if pname in expnames:
            keep.append(list(expnames).index(pname))
            pnamekeep.append(p)
    sequences = np.array(sequences)[keep]
    kweight = kweight[pnamekeep]
    pnames = pnames[pnamekeep]
    


    
outadd = ''    
if '--in-kmersignificance' in sys.argv:
    import scipy.stats as sct
    outadd = '_in-ksig'
    print 'in-kmersignificance', np.shape(kweight)
    def pval(kcount, mean, std,jj):
        print jj
        pvalsout = np.zeros(len(kcount))
        for y, kjc in enumerate(kcount):
            pvalsout[y] = -np.log10(sct.norm.sf(kjc,mean[y],std[y]))
        return jj, pvalsout
    
    kmean = np.mean(Ppred, axis = 0)
    kstd = np.std(Ppred, axis = 0)
    print kmean, kstd
    results = Parallel(n_jobs=num_cores)(delayed(pval)(kweight[j], kmean, kstd, j) for j in range(len(kweight)))
    print results
    for i in  range(len(kweight)):
        kweight[results[i][0]] = results[i][1]


elif '--between-kmersignificance' in sys.argv:
    import scipy.stats as sct
    outadd = '_btw-ksig'
    print 'between-kmersignificance'
    kmean = np.mean(Ppred.flatten())
    kstd = np.std(Ppred.flatten())
    print kmean, kstd

    def pval(kcount, mean, std, jj):
        print jj
        pvalsout = np.zeros(len(kcount))
        for y, kjc in enumerate(kcount):
            sfval = sct.norm.sf(kjc,mean,std)
            sfval = max(1e-300, sfval)
            pvalsout[y] = -np.log10(sfval)
        return jj, pvalsout
    
    results = Parallel(n_jobs=num_cores)(delayed(pval)(kweight[j], kmean, kstd, j) for j in range(len(kweight)))
    print results
    for i in  range(len(kweight)):
        kweight[results[i][0]] = results[i][1]




if '--cutoffweight' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffweight')+1])
    outadd += '-cut'+str(wkcut)
    kweight[kweight <= wkcut] = 0
    print 'Features left', np.sum(kweight > 0, axis = 1)

elif '--cutoffmean' in sys.argv:
    outadd += '-cut-mean'
    kweight[kweight <= np.mean(Ppred)] = 0
    print 'Features left', np.sum(kweight > 0, axis = 1)

elif '--cutofftop' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutofftop')+1])
    outadd += '-cuttop'+str(wkcut)
    wkcut = np.sort(Ppred.flatten())[-int(wkcut*float(len(Ppred.flatten())))]
    kweight[kweight <= wkcut] = 0
    
elif '--cutoffcentre' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffcentre')+1])
    outadd += '-cutcentre'+str(wkcut)
    wkcut = np.sort(Ppred.flatten())[-int(wkcut*float(len(Ppred.flatten())))]
    wlcut = np.sort(Ppred.flatten())[int(wkcut*float(len(Ppred.flatten())))]
    kweight[(kweight <= wkcut) & (kweight >= wlcut)] = 0

elif '--cutoffsignificanttop' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffsignificanttop')+1])
    outadd += '-cutsigniftop'+str(wkcut)
    meanall = np.mean(Ppred.flatten())
    stdall = np.std(Ppred.flatten())
    kweight[kweight <= wkcut * std + meanall] = 0
    
elif '--cutoffsignificantcentre' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffsignificantcentre')+1])
    outadd += '-cutsignifcentre'+str(wkcut)
    meanall = np.mean(Ppred.flatten())
    stdall = np.std(Ppred.flatten())
    kweight[(kweight <= wkcut * std + meanall) & (kweight >= -wkcut * std + meanall)] = 0
    


elif '--cutofftop-each' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutofftop-each')+1])
    cutlen = int(wkcut * float(len(Ppred[0])))
    print 'Keep top', cutlen
    outadd += '-cuttop-each'+str(wkcut)
    for kw, kwpfeat in enumerate(kweight):
        kwpfeat[kwpfeat <= np.sort(kwpfeat)[-cutlen]] = 0
        kweight[kw] = kwpfeat
        
elif '--cutoffcentre-each' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffcentre-each')+1])
    cutlen = int(wkcut * float(len(Ppred[0])))
    print 'Keep extremes', cutlen
    outadd += '-cutcetnre-each'+str(wkcut)
    for kw, kwpfeat in enumerate(kweight):
        kwpfeat[(kwpfeat <= np.sort(kwpfeat)[-cutlen]) & (kwpfeat >= np.sort(kwpfeat)[cutlen])] = 0
        kweight[kw] = kwpfeat
        
elif '--cutoffmeantop-each' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffmeantop-each')+1])
    outadd += '-cuttopmean-each'+str(wkcut)
    for kw, kwpfeat in enumerate(kweight):
        checklen = kwpfeat[(kwpfeat - np.mean(kwpfeat)) > 0]
        cutlen = int(wkcut * float(len(checklen)))
        print 'Keep top', cutlen
        kwpfeat[kwpfeat <= np.sort(kwpfeat)[-cutlen]] = 0
        kweight[kw] = kwpfeat

elif '--cutoffsignificanttop-each' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffsignificanttop-each')+1])
    outadd += '-cuttopsig-each'+str(wkcut)
    for kw, kwpfeat in enumerate(kweight):
        meankw = np.mean(kwpfeat)
        stdkw = np.std(kwpfeat)
        kwpfeat[kwpfeat <= meankw + wkcut * stdkw] = 0
        print 'Kept', len(np.nonzero(kwpfeat)[0]), meankw, wkcut*stdkw, np.amin(kwpfeat), np.amax(kwpfeat)
        kweight[kw] = kwpfeat
        
elif '--cutoffsignificantcentre-each' in sys.argv:
    wkcut = float(sys.argv[sys.argv.index('--cutoffsignificantcentre-each')+1])
    outadd += '-cutsigcentre-each'+str(wkcut)
    for kw, kwpfeat in enumerate(kweight):
        meankw = np.mean(kwpfeat)
        stdkw = np.std(kwpfeat)
        kwpfeat[(kwpfeat <= meankw + wkcut * stdkw) & (kwpfeat >= meankw - wkcut * stdkw)] = 0
        print 'Kept', pnames[kw], len(np.nonzero(kwpfeat)[0]), meankw, wkcut*stdkw, np.amin(kwpfeat), np.amax(kwpfeat)
        kweight[kw] = kwpfeat


if '--original_kmerfeatures' in sys.argv:
    kweight = kweight[:, preckeep]
    for p, pname in enumerate(pnames):
        op = list(ornames).index(pname)
        kweight[p] = kweight[p] * (Porig[:, op] > 0).astype(float)
        print len(np.nonzero(kweight[p])[0]), np.amin(kweight[p]), np.amax(kweight[p])


if '--multiprocessing' in sys.argv:
    num_cores = int(sys.argv[sys.argv.index('--multiprocessing')+1])
    results = Parallel(n_jobs=num_cores)(delayed(map2seq)(kweight[t], Pfeatures, sequences[t], t) for t in range(0,len(pnames)))

    resultsindex = [results[i][0] for i in range(len(results))]

    sortsort = np.argsort(resultsindex)
    sequenceweights = []
    for so in sortsort:
        sequenceweights.append(results[so][1])
else:
    sequenceweights =[]
    for t in range(len(pnames)):
        it, seqweight = map2seq(kweight[t], Pfeatures, sequences[t], t)
        sequenceweights.append(seqweight)
   



if '--savesequencemap' in sys.argv:
    soutname = sys.argv[sys.argv.index('--savesequencemap')+1]
    print "Saved sequence maps as", soutname
    sobj = open(os.path.splitext(soutname)[0]+'.map', 'w+')
    for t, testprot in enumerate(pnames):
        sobj.write('%15s' % testprot)
    sobj.write('\n')
    for i in range(len(max(sequences,key=len))):
        for t, testprot in enumerate(pnames):
            if i < len(sequences[t]):
                sobj.write('%15s' % (str(i)+' '+sequences[t][i]+' '+str(round(sequenceweights[t][i],3))))
            else:
                sobj.write('%15s' % '')
        sobj.write('\n')
else:
    soutname = startout+os.path.splitext(os.path.split(common_name)[1])[0]+'_kmer2sequence'+outadd+'.fasta'
    print "Saved sequence maps as", soutname
    sobj = open(soutname, 'w')
    for t, testprot in enumerate(pnames):
        #print '\n'+','.join(np.array(sequenceweights[t], dtype = str))
        sobj.write('>'+testprot+'\n'+sequences[t]+'\n'+','.join(np.around(sequenceweights[t],3).astype(str))+'\n')
    

