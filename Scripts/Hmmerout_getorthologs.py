import numpy as np
import sys, os 
from joblib import Parallel, delayed
import multiprocessing

def readinhmmer(hmfile):
    obj = open(hmfile, 'r').readlines()
    
    queries = []
    fullsequences = [] # Sequence with gaps in templates sequence from ortholog domain
    fulldomainsequences = [] # Concatenated sequence or ortholog protein
    sequences = [] # Sequences aligned to the input sequence
    seqscores = []
    seqids = []
    seqnames = []
    domainnames = []
    seqaccrons = []
    species = []
    
    # Sets for each template sequence
    tfullsequences = [] # Sequence with gaps in templates sequence from ortholog domain
    tfulldomainsequences = [] # Concatenated sequence or ortholog protein
    tsequences = [] # Sequences aligned to the input sequence
    tseqscores = []
    tseqids = []
    tseqnames = []
    tseqaccrons = []
    tspecies = []
    tdomainnames = []
    fsequences = []
    newprot = False
    read_out = False
    protnum = -1
    domnum = -1
    for l, line in enumerate(obj):
        # read in Query protein
        if line[:6] == 'Query:' and obj[l+6].strip()[1:8] != 'No hits':
            print line
            origname = line.split()[1]
            origlen = int(line.strip().split()[2][3:-1])
            fsequences = []
        elif line[:6] == 'Query:' and obj[l+6].strip()[1:8] == 'No hits':
            print line
            print obj[l+6]
            fsequences = []
            
        
        # Read out info from sequences
        if read_out:
                    
            if line[:2] == '>>' and obj[l+1].strip()[:14]!='[No individual':
                #print line
                protnum += 1
                if protnum > 0:
                    #print protnum, len(fsequences)
                    if domnum == 0:
                        fsequences = [domseq.replace('X', '-')]
                    tfullsequences.append('*'.join(np.array(fsequences)))
                    tspecies.append(orthspec)
                    tseqaccrons.append(orthoacc)
                    tseqnames.append(orthoname)
                # read in name of aligned sequence
                orthoacc = line.split()[1]
                if 'OS=' in line:
                    orthoname = line.split(orthoacc)[1].split('OS=')[0]
                    orthspec = line.split(orthoacc)[1].split('OS=')[1].strip()
                elif 'sp|' in line:
                    orthoname = line.split('|')[1]
                    orthspec = line.split('|')[2]
                elif '||' in line:
                    orthoname = line.split('||')[1].split('|')[1].strip()
                    orthspec = line.split('||')[1].split('|')[0].strip()
                elif '|' in line:
                    orthoname = line.split('|')[1]
                    orthspec = line.split('|')[2].split('__')[0]
                else:
                    orthoname = line.strip()
                    orthspec = line.strip()
                    #print orthspec, orthoname
                # define whole sequence variables that need to be completed
                startprot = 0
                fsequences = []
            
            if line[:11] == '  == domain':
                domnum += 1
                if domnum > 0:
                    #print 'Dnum', domnum
                    tsequences.append(''.join(alignseq).replace('X', '-'))
                    tfulldomainsequences.append(domseq.replace('X', '').replace('-', ''))
                    tseqscores.append(bitscore)
                    tseqids.append(evals)
                    tdomainnames.append(domname)
                    fsequences.append(domseq.replace('X', '').replace('-', ''))
                    
                evals = float(line.split()[-1])
                bitscore = float(line.split()[4])
                domname = orthoacc+'_d'+line.split()[2]
                alignseq = np.array(['-' for o in range(origlen)])
                domseq = ''
            
            if not line.strip():
                readortho = False
                continue
            else:
                #print line.split()[0], origname, orthoacc
                if line.split()[0] == origname:
                    #print line
                    readortho = True
                    line = line.split()
                    location = [int(line[1])-1, int(line[3])]
                    #print location, len(line[2])
                    if '.' in line[2]:
                        masking = True
                        mask = np.array(list(line[2])) != '.'
                        if np.sum(mask) != location[1] - location[0]:
                            location[1] = location[0] + int(np.sum(mask))
                    else:
                        masking = False
                elif line.split()[0] == orthoacc and readortho:
                    #print line
                    if masking:
                        lineseq = np.array(list(line.split()[2]))
                        alignseq[location[0]:location[1]] = lineseq[mask]
                    else:
                        #print len(line.split()[2])
                        alignseq[location[0]:location[1]] = list(line.split()[2])
                    domseq += line.split()[2]
                    readortho = False
        if line[:36] == 'Internal pipeline statistics summary' and len(fsequences) > 0:
            tsequences.append(''.join(alignseq).replace('X', '-'))
            tfulldomainsequences.append(domseq.replace('X', '').replace('-', ''))
            fsequences.append(domseq.replace('X', '').replace('-', ''))
            tseqscores.append(bitscore)
            tseqids.append(evals)
            tdomainnames.append(domname)                    
            tfullsequences.append('*'.join(np.array(fsequences)))
            tspecies.append(orthspec)
            tseqaccrons.append(orthoacc)
            tseqnames.append(orthoname)
            newprot = True
            print len(tfulldomainsequences), len(tfullsequences), len(tsequences), len(tspecies), len(tseqaccrons), len(tseqnames), len(tdomainnames)
        elif line[:36] == 'Internal pipeline statistics summary' and len(fsequences) == 0:
            newprot = True
            print l, line
        # Determine start of sequences in file
        if line[:17] == 'Domain annotation' and obj[l+1][:2] == '>>':
            read_out = True
            startprot = -1
        if newprot and len(tfullsequences) > 0:
            queries.append(origname)
            fullsequences.append(tfullsequences)
            fulldomainsequences.append(tfulldomainsequences)
            sequences.append(tsequences)
            seqscores.append(tseqscores)
            seqids.append(tseqids)
            seqnames.append(tseqnames)
            seqaccrons.append(tseqaccrons)
            species.append(tspecies)
            domainnames.append(tdomainnames)
            
            tfullsequences = [] # Sequence with gaps in templates sequence from ortholog domain
            tfulldomainsequences = [] # Concatenated sequence or ortholog protein
            tsequences = [] # Sequences aligned to the input sequence
            tseqscores = []
            tseqids = []
            tseqnames = []
            tseqaccrons = []
            tspecies = []
            tdomainnames = []
            read_out = False
            newprot = False
            protnum = -1
            domnum = -1
            startprot = 0
            fsequences = []
        elif newprot and len(fsequences) == 0:
            print 'Skipped'
            tfullsequences = [] # Sequence with gaps in templates sequence from ortholog domain
            tfulldomainsequences = [] # Concatenated sequence or ortholog protein
            tsequences = [] # Sequences aligned to the input sequence
            tseqscores = []
            tseqids = []
            tseqnames = []
            tseqaccrons = []
            tspecies = []
            tdomainnames = []
            read_out = False
            newprot = False
            protnum = -1
            domnum = -1
            startprot = 0
            fsequences = []
    
    return  queries, sequences, fullsequences, fulldomainsequences, seqscores, seqids, seqnames, domainnames, seqaccrons, species
                
def getidmsa(smsa):
    msaid = np.zeros((len(smsa), len(smsa)))
    lenmsamat = np.ones((len(smsa), len(smsa)))
    for s in range(len(smsa)):
        mulsa = smsa[s]
        for t in range(s, len(smsa)):
            tulsa = smsa[t]
            lenmsa = 0.
            idmsa = 0.
            for p in range(len(mulsa)):
                pos = mulsa[p]
                if pos != '-' or tulsa[p] != '-':
                    lenmsa += 1.
                    if pos == tulsa[p]:
                        idmsa +=1.
            msaid[s,t] = msaid[t,s] = idmsa
            lenmsamat[s,t] = lenmsamat[t,s] = lenmsa
    return msaid/ lenmsamat

def getdomainidmsa(smsa, domregion):
    msaid = np.zeros((len(smsa), len(smsa)))
    msadomid = np.zeros((len(domregion), len(smsa), len(smsa)))
    for s in range(len(smsa)):
        mulsa = smsa[s]
        for t in range(s, len(smsa)):
            tulsa = smsa[t]
            lenmsa = 0.
            idmsa = 0.
            for d, dregion in enumerate(domregion):
                iddom = 0.
                for p in range(dregion[0], dregion[1]):
                    posA = mulsa[p]
                    posB = tulsa[p]
                    if posA != '-' or posB != '-':
                        lenmsa += 1.
                        if posA == posB:
                            idmsa += 1.
                            iddom += 1.
                msadomid[d,s,t] = msadomid[d, t, s] = iddom/float(dregion[1]-dregion[0])
            msaid[s,t] = msaid[t,s] = float(idmsa)/float(lenmsa) 
    return msaid, msadomid

def threshclust(distcut, clstnames, cteringmat, cctype):
        clustpos = np.arange(len(clstnames),dtype = int)
        db = -np.ones(len(clstnames), dtype = int)
        clustcentre = []
        clusters = []
        cluspscores = []
        c = 0
        while True:
            setlens, closesets = clustersize(cteringmat, distcut)
            lcl = np.argmax(setlens)
            clustcentre.append(clustpos[lcl])
            aclust = addtoclust(closesets[lcl], cteringmat, lcl, distcut, cctype)
            clusters.append(clustpos[aclust])
            eescores = cteringmat[aclust]
            cluspscores.append(np.sum(eescores[:, aclust], axis = 0))
            db[clustpos[aclust]] = c
            c+=1
            clustpos = np.delete(clustpos, aclust)
            cteringmat = np.delete(cteringmat, aclust, axis = 0)
            cteringmat = np.delete(cteringmat, aclust, axis = 1)
            if len(cteringmat) == 0:
                break
        return clustcentre, clusters, db, cluspscores
    
def clustersize(cmat, dcut):
    closesets = []
    setlens = []
    for i in range(len(cmat)):
        closesets.append(np.where(cmat[i] >= dcut)[0])
        setlens.append(len(closesets[-1]))
    setlens = np.array(setlens)
    return setlens, closesets
    
def addtoclust(cset, cmat, lc, dcut, wtype):
    if wtype == 'centre':
        aset = cset
    elif wtype == 'weak':
        dsetsize = 1
        while dsetsize > 0:
            osetsize = len(cset)
            cnsets = np.where(cmat[cset] > dcut)[1]
            addset = []
            cset = np.unique(np.append(cset, cnsets))
            if len(cset) > 1:
                highidnp = np.argsort(cmat[cset], axis = 1)[:,-2]
                for cne, cnentry in enumerate(cset):
                    if highidnp[cne] in cset:
                        addset.append(cne)
                cset = cset[addset]
            nsetsize = len(cset)
            dsetsize = nsetsize - osetsize 
            #print cset, nsetsize, osetsize, dsetsize
        aset = cset    
    elif wtype == 'weak-overlap':
        # Determines the overlap with the current and the other clusters and sorts by this
        dsetsize = 1
        while dsetsize > 0:
            osetsize = len(cset)
            cnsets = np.where(cmat[cset] > dcut)[1]
            addset = []
            cset = np.unique(np.append(cset, cnsets))
            if len(cset) > 1:
                # define overlap with different clusters
                for cne, cnentry in enumerate(cset):
                    potentialcluster = np.where(cmat[cnentry] > dcut)[0]
                    potentclustercomp = np.where(cmat[potentialcluster] > dcut)
                    potentialoverlap = []
                    for po in range(len(potentialcluster)):
                        potento = len(np.intersect1d(potentialcluster, potentclustercomp[1][potentclustercomp[0] == po]))
                        potentialoverlap.append(potento)
                    potentialoverlap = np.array(potentialoverlap)
                    largest_over = np.amax(potentialoverlap)
                    if lc in potentialcluster[potentialoverlap == largest_over]:
                        addset.append(cne)
                cset = cset[addset]
            nsetsize = len(cset)
            dsetsize = nsetsize - osetsize
        aset = cset  
        
    elif wtype == 'strong':
        aset = [lc]
        csort = np.argsort(-cmat[lc,cset])
        cset = cset[csort]
        for cs in cset:
            addto = True
            if cs == lc:
                addto = False
            else:
                for closen in cmat[cs,aset]:
                    if closen < dcut:
                        addto = False
            if addto:
                aset.append(cs)
    return aset
        
def getlenalign(alignseq):
    lenalign = []
    for alseq in alignseq:
        lenalign.append(len(np.where(np.array(list(alseq)) != '-')[0]))
    return np.array(lenalign)

def getdomaincoverage(domainpos, alignseq):
    coverage = []
    for alseq in alignseq:
        cov = []
        for domainp in domainpos:
            lenaldom = domainp[1] - domainp[0]
            lencov = float(len(np.where(np.array(list(alseq))[domainp[0] :domainp[1]] != '-')[0]))
            cov.append(lencov/float(lenaldom))
        coverage.append(cov)
    return np.array(coverage)

def removeseqs(idmat, maxsim, minhom):
    keep = np.where(idmat[0] >= minhom)[0]
    rm = []
    for ki, k in enumerate(keep):
        for li in range(ki + 1, len(keep)):
            l = keep[li]
            if idmat[k,l] > maxsim:
                rm.append(l)
    outkeep = np.setdiff1d(keep, rm)
    return outkeep

def removeseqsgetid(seqtest, maxsim, minhom, numcores):
    outkeep = np.arange(len(seqtest))
    print 'Start', len(outkeep)
    unseqs, unseqind = np.unique(seqtest, return_index= True)
    print 'Unique', len(unseqs)
    remo = np.delete(np.arange(len(seqtest)), unseqind)
    seqtest = np.delete(seqtest, remo)
    outkeep = np.delete(outkeep, remo)
    print 'Unique in left', len(outkeep)
    idvec = np.ones(len(seqtest))
    
    lennorm = float(len(seqtest[0]))
    oseq = seqtest[0]
    remo = []
    for s, seq in enumerate(seqtest[1:]):
        # only adjust lennorm when identity should only be measured in part without gaps -, ends up with high identity fragments
            #lennorm = lenfunc(seq)
        cidf, s = idseq(oseq, seq, s)
        cid = cidf/lennorm
        #print cid
        idvec[s+1] = cid
        if cid > maxsim:
            remo += [s + 1]
        elif cid < minhom:
            remo += [s + 1]
    idvec = np.delete(idvec, remo)
    seqtest = np.delete(seqtest, remo)
    outkeep = np.delete(outkeep, remo)
    print 'Left', len(outkeep)
    s = 0
    
    idmat = np.ones((len(seqtest), len(seqtest)))
    idmat[0] = idvec
    idmat[:, 0] = idvec
    #lennorm = float(len(seqtest[0]))
    while True:
        remo = []
        s += 1
        if s% 100 == 0:
            print s, len(seqtest)
        if s == len(seqtest):
            break
        seq = seqtest[s]
        lennorm = len(seq)
        
        #if numcores == 1:
        for t in range(s+1, len(seqtest)):
                seqb = seqtest[t]
                cidf, t = idseq(seq, seqb,t)
                cid = cidf/lennorm
                idmat[s,t] = idmat[t, s] = cid
                if cid > maxsim:
                    remo += [t]    
        #else:
            #results = Parallel(n_jobs=numcores)(delayed(idseq)(seq, seqtest[t], t) for t in range(s+1,len(seqtest)))
            #for cidf, t in results:
                #cid = cidf/lennorm
                #idmat[s,t] = idmat[t, s] = cid
                #if cid > maxsim:
                    #remo += [t]
        idmat = np.delete(idmat, remo, axis = 0)
        idmat = np.delete(idmat, remo, axis = 1)
        seqtest = np.delete(seqtest, remo)
        outkeep = np.delete(outkeep, remo)
    print 'Final', len(outkeep)
    return seqtest, idmat, outkeep

#@jit(nopython=True, parallel=True)            
def idseq(mulsa, tulsa, t):
    idmsa = 0.
    for p in range(len(mulsa)):
        pos = mulsa[p]
        if pos == tulsa[p]:
            idmsa +=1.
    return idmsa, t

def lenfunc(someseq):
    gaploc = np.where(np.array(list(someseq)) != '-')[0]
    return float(gaploc[-1]-gaploc[0])
            
if __name__ == '__main__':
    if '--multiprocessing' in sys.argv:
        numcores = int(sys.argv[sys.argv.index('--multiprocessing')+1])
    else:
        numcores = 1
    
    hmmerfile = sys.argv[1]

    queries, sequences, fullsequences, fulldomainsequences, seqscores, seqids, seqnames, domainnames, seqaccrons, species = readinhmmer(hmmerfile)
    
    #print queries, len(queries)

    fastafile = sys.argv[2]
    faslines = open(fastafile, 'r').readlines()
    fastaseq = []
    fastaname = []
    for f, fline in enumerate(faslines):
        if fline[0] == '>':
            fastaname.append(fline.strip()[1:])
            fastaseq.append(faslines[f+1].strip().replace('X', '-'))
    
    maxsim = float(sys.argv[3])
    minhom = float(sys.argv[4])
    outfolder = sys.argv[5]
    print 'readin'
    if os.path.splitext(outfolder)[-1] == '.fasta':
        print 'Single file out:', outfolder
        singlefile = True
        mobj = open(os.path.splitext(outfolder)[0]+'_homologs_maxsim'+str(maxsim)+'_minhom'+str(minhom)+'.fasta', 'w')
    else:
        singlefile = False
    # First calculate only identy to mastersequence and delete all sequences that are too far
    # unite identity finding with sequence selection to speed up

    for f, qfasta in enumerate(fastaname):
        print f, qfasta
        if qfasta in queries:
            q = queries.index(qfasta)
            query = queries[q]
            qfastseq = fastaseq[f]
            qseqs = np.append([qfastseq], sequences[q])
            if '||' in query:
                query = query.split('||')[0]
            elif '|' in query:
                query = query.rsplit('|',1)[0].replace('|','_')+"__"+query.rsplit("__",1)[1]
            
            #print qseqs, maxsim, minhom
            final_seq, idents, keepset = removeseqsgetid(qseqs, maxsim, minhom, numcores)
            print 'Kept sequences', len(keepset)
            if len(idents) > 1:
                print 'Min id to target', np.amin(idents[0]), 'max id to target', np.amax(idents[0][1:])
            print 'Min between', np.amin(idents), 'max between', np.amax(np.triu(idents, k=1))
            if singlefile:
                #mobj.write('>'+query+'\n'+qfastseq+'\n')
                if len(keepset[1:]) > 0:
                    for kkep in keepset[1:]:
                        mobj.write('>'+qfasta+'||_'+domainnames[q][kkep-1]+'\n'+sequences[q][kkep-1]+'\n')
            else:
                msaobj = open(outfolder+'/'+query+'_homologs_maxsim'+str(maxsim)+'_minhom'+str(minhom)+'_msa.fasta', 'w')
                fasobj = open(outfolder+'/'+query+'_homologs_maxsim'+str(maxsim)+'_minhom'+str(minhom)+'.fasta', 'w')
                msaobj.write('>'+qfasta+'\n'+qfastseq+'\n')
                fasobj.write('>'+qfasta+'\n'+qfastseq+'\n')
                if len(keepset[1:]) > 0:
                    for kkep in keepset[1:]:
                        msaobj.write('>'+domainnames[q][kkep-1]+'\n'+sequences[q][kkep-1]+'\n')
                        fasobj.write('>'+domainnames[q][kkep-1]+'\n'+fulldomainsequences[q][kkep-1]+'\n')
                msaobj.close()
                fasobj.close()
                if '--saveidmat' in sys.argv:
                    np.savetxt(outfolder+'/'+query+'_homologs_maxsim'+str(maxsim)+'_minhom'+str(minhom)+'.idmat', idents, fmt = '%1.3f', header=query+' '+' '.join(np.array(domainnames[q])[keepset[1:]]))
        else:
            print 'Not included'
            
            
        




                
