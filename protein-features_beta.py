#protein-features.py
import numpy as np
import sys, os
from itertools import imap
from operator import eq
from joblib import Parallel, delayed
import multiprocessing
if '--mprocessing' in sys.argv:
    mprocessing = True
    num_cores = int(sys.argv[sys.argv.index('--mprocessing')+1])
else:
    mprocessing = False
    num_cores = 1

def dposloc(dpstring):
    for dp in range(len(dpstring)):
        if not dpstring[-1-dp].isdigit():
            break
    return -dp


def blosumsub(kmerfeat, simsmomat, maxdiff,kmeroutname):
    # corrects for a frequency of mutations in potential homolgs
    # f.e. 10% of AA be mutated == 90 sequence identity
    if '--mutation_frequency' in sys.argv:
        mfreq = float(sys.argv[sys.argv.index('--mutation_frequency')+1])
        kmeroutname += '_mutfreq'+str(mfreq)
        mfreqd = 1.-mfreq
    else:
        # mutation frequency is not respected ==>
        mfreq = 1.
        mfreqd = 1.
    
    # only allow mutations that have a higher blosum score than:
    if '--minblos' in sys.argv:
        minblos = float(sys.argv[sys.argv.index('--minblos')+1])
        kmeroutname += '_minblosum'+str(minblos)
    else:
        minblos = -20
    
    def blosim(ksim1, ksim2):
        simscore = []
        aextype = []
        for aap in zip(ksim1, ksim2):
            respo1 = residues.index(aap[0])
            respo2 = residues.index(aap[1])
            bloss = simsmomat[respo1,respo2]
            #print aap, respo1, respo2, bloss
            if bloss < minblos:
                simscore = []
                return False, aextype, simscore
            else:
                aextype.append(respo2)
                simscore.append(bloss)
        return True, aextype, simscore
    
    def eq(ksim1, ksim2):
        simscore = []
        for aap in zip(ksim1, ksim2):
            simscore.append(aap[0] == aap[1])
        return simscore
    
    if '--gausnorm' in sys.argv:
        kmeroutname += '_gausnorm'
        def normblos(bvals, bvalsmax, bpos, mpos):
            bvalsmax = np.sum(bvalsmax)
            bvals = np.sum(bvals)
            bvals = np.exp(-(bvals-bvalsmax)**2/bvalsmax)
            return bvals
    elif '--quadraticnorm' in sys.argv:
        kmeroutname += '_quadnorm'
        def normblos(bvals, bvalsmax, bpos, mpos):
            bvalsmax = np.sum(bvalsmax)**2
            bvals = np.sum(bvals)**2
            bvals = bvals/bvalsmax
            return bvals
    elif '--lognorm' in sys.argv:
        kmeroutname += '_lognorm'
        nummut = simsmomat[:-1,:-1] >= minblos
        nummut = nummut * backfreq
        nummut = nummut/np.sum(nummut, axis = 1)
        if '--onlymfreq' in sys.argv:
            kmeroutname += '-onlymfreq'
            #print nummut
            def qij(q,w, y):
                return 1.*nummut[y,q]
            def qii(q,w, y):
                return 1.*nummut[y,q]
        else:
            def qij(q,w, y):
                return  2.*2.**(w/2.) * nummut[y,q]
            def qii(q,w, y):
                return  2.**(w/2.) * nummut[y,q]
        
        def normblos(bvals, bvalsmax, bpos, mpos):
            bscore = 1.
            for fg, [f,g] in enumerate(zip(bvals, bvalsmax)):
                if f == g:
                    bscore *= mfreqd * qii(bpos[fg],g, mpos[fg])
                else:
                    bscore *= mfreq * qij(bpos[fg],f,mpos[fg])
                #print bscore
            return bscore
    elif '--nonorm' in sys.argv:
        kmeroutname += '_nonorm'
        def normblos(bvals, bvalsmax, bpos, mpos):
            return 1.
    else:
        kmeroutname += '_linearnorm'
        def normblos(bvals, bvalsmax, bpos, mpos):
            bvalsmax = np.sum(bvalsmax)
            bvals = np.sum(bvals)
            bvals = bvals/bvalsmax
            return bvals                
    
    # define the weights given to each kmer and its variants. 
    print 'Generate kmer substitution matrix...'
    subpos = []
    subval = []
    #kmersfeat contain all sets of kmers, 1mer, 2mer, 3mer, ...
    if len(np.shape(kmerfeat)) == 1:
        kmerfeat = np.array([kmerfeat])
    for ka, kmers in enumerate(kmerfeat):
        if ka > 0:
            kaddlen = len(kmerfeat[ka-1])
        else:
            kaddlen = 0
        if len(kmers[0]) - maxdiff >= 2:
            lensim = len(kmers[0]) - maxdiff
            # for each kmer in the set determine similar kmers, their position and their weight
            for kmer in kmers:
                # get score for alignment of right kmer
                det, aapos, kmerscore = blosim(kmer, kmer)
                kpos = []
                kvals = []
                # go through all other kmers and compare
                for ks, simkmer in enumerate(kmers):
                    # check if minimum positions are
                    if sum(eq( kmer, simkmer)) >= lensim:
                        addreal, simaapos, simkweight = blosim(kmer, simkmer)
                        if addreal:
                            # norm the weights
                            nweight = normblos(simkweight, kmerscore, simaapos, aapos)
                            #print addreal, kmer, simkmer, simaapos, simkweight
                            #print nweight
                            kpos.append(ks + kaddlen)
                            kvals.append(nweight)
                subpos.append(kpos)
                subval.append(kvals)
        else:
            for i in range(kaddlen, kaddlen+len(kmers)):
                subpos.append([i])
                subval.append([1])
    return subpos, subval, kmeroutname


if '--fastafmt' in sys.argv:
    sequencefile = sys.argv[sys.argv.index("--fastafmt")+1]
    sobj = open(sequencefile, 'r').readlines()
    protnames = []
    expnames = []
    sequence = []
    score = []
    scoretype = '1'
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
                #if 'X' in objseq:
                    #objseq = objseq.replace("X", "*")
                if '*' in objseq:
                    objseq = objseq.split('*')
                else:
                    objseq = [objseq]
                sequence.append(objseq)
                sco = []
                for seq in sequence[-1]:
                    sco.append(np.ones(len(seq)))
                score.append(sco)

elif '--domainfasta' in sys.argv:
    fi = sys.argv[sys.argv.index("--domainfasta")+1]
    fo = open(fi, 'r').readlines()
    names = []
    seqs = []
    dnums = []
    for l, line in enumerate(fo):
        if line[0] == '>':
            names.append(line.split('__')[0])
            seqs.append(fo[l+1].strip())
            dnums.append(int(line.split('_')[-1]))
    names = np.array(names)
    seqs = np.array(seqs)
    unames = np.unique(names)
    dnums = np.array(dnums)
    protnames = []
    expnames = []
    sequence = []
    score = []
    scoretype = '1'
    for u in range(len(unames)):
        uname = unames[u]
        # clean uname if necessary
        if '||' in uname:
            protnames.append(uname.split("||")[0][1:])
            expnames.append(uname.split("||")[1])
        else:
            protnames.append(uname[1:])
            expnames.append(uname[1:])
        
        objnums = dnums[names == unames[u]]
        objseq = ['' for ju in range(np.amax(objnums)+1)]
        cobjseq = seqs[names == unames[u]]
        for o, obn in enumerate(objnums):
            objseq[obn] = cobjseq[o]
        #print objseq
        sequence.append(objseq)
        
        sco = []
        for seq in sequence[-1]:
            sco.append(np.ones(len(seq)))
        score.append(sco)
    #sys.exit()
    #print sequence[213]
    #print expnames
        
                
elif "--infofile" in sys.argv:
    ### .info files
    sequencefile = sys.argv[sys.argv.index("--infofile")+1]
    sobj = open(sequencefile, 'r')
    protnames = []
    expnames = []
    sequence = []
    score = []
    scoretype = '1'
    for i, sline in enumerate(sobj):
            sline = sline.strip().split()
            if '||' in sline[0]:
                if sline[0][0] == '>':
                    protnames.append(sline[0].split("||")[0][1:])
                else:
                    protnames.append(sline[0].split("||")[0])
                expnames.append(sline[0].split("||")[1])
            else:
                protnames.append(sline[0])
                expnames.append(sline[0])
            ndom = int(sline[1])
            sequence.append(sline[-ndom:])
            sco = []
            for seq in sequence[-1]:
                sco.append(np.ones(len(seq)))
            score.append(sco)
    if '_' in protnames[0]:
        for pm, protname in enumerate(protnames):
            protnames[pm] =protname.split('_')[0]
    


elif '--sequencescores' in sys.argv:
    sequencefile = sys.argv[sys.argv.index("--sequencescores")+1]
    scoretype = sys.argv[sys.argv.index("--sequencescores")+2]
    sobj = open(sequencefile, 'r')
    protnames = []
    expnames = []
    sequence = []
    score = []
    for i, sline in enumerate(sobj):
            if len(sline.split()) == 2:
                sline = sline.strip().split()
                protnames.append(sline[0].split("||")[0])
                expnames.append(sline[0].split("||")[1])
                ndom = int(len(sline[-1].split('*')))
                sequence.append(sline[-1].split('*'))
            else:
                sco = sline.strip().split(',*,')
                scor = []
                for sc in sco:
                    scor.append(np.array(sc.split(','), dtype = float))
                score.append(scor)
    if scoretype == 'binary':
        scoretype = 'max'
        cut = float(sys.argv[sys.argv.index("--sequencescores")+3])
        for s, scor in enumerate(score):
            for t, tcor in enumerate(scor):
                score[s][t][tcor >= cut] = 1.
                score[s][t][tcor < cut] = 0.
            
    
else:
    print "Please provide one of the possible file formats"
    sys.exit()

def outscore(Prob, typ):
    if typ == 'max':
        pout = np.amax(Prob)
    elif typ == 'mean':
        pout = np.mean(Prob)
    elif typ == 'sum':
        pout == np.sum(Prob)
    elif typ == '1':
        pout = 1
    return pout


if "--k" in sys.argv:
    k = int(sys.argv[sys.argv.index("--k")+1])
else: 
    print "Please give the length of the kmer" 
    sys.exit()
    
allk = False
if "--allk" in sys.argv:
    allk = True
    start = 1
else:
    start = k

if '--startk' in sys.argv:
    allk = True
    start = int(sys.argv[sys.argv.index('--startk')+1])

if "--outname" in sys.argv:
    outname = sys.argv[sys.argv.index("--outname")+1]
else:
    print 'Please define outname'
    sys.exit()

if '--orthologs' in sys.argv:
    # Proteinorthologdomain-gapped.info
    orthseqs = sys.argv[sys.argv.index("--orthologs")+1]
    orthofmt = os.path.splitext(orthseqs)[1]
    # before one  had (sequence, domains), now it's going to be (seq, dom, orthologs)
    oobj = open(orthseqs, 'r').readlines()
    onames = []
    odomseqs = []
    if orthofmt == '.info':
        for line in oobj:
            line = line.strip().split()
            onumdom = int(line[1])
            odomseqs.append(line[-onumdom:])
            onames.append(line[0].split("_",1)[0].split("||")[1])
        # add ortholg sequence to sequence array
        onames = np.array(onames)
        for e, ename in enumerate(expnames):
            if ename in onames:
                sequence[e] = [sequence[e]]
                for t in np.where(onames == ename)[0]:
                    sequence[e].append(odomseqs[t])
                sequence[e] = np.array(sequence[e]).T
            else:
                sequence[e] = np.array([sequence[e]]).T
    if orthofmt == '.fasta':
        for l, line in enumerate(oobj):
            if line[0] == '>':
                odomseqs.append(oobj[l+1].strip().split('*'))
                horigname = line[1:].split("||_",1)[0]
                if '||' in horigname:
                    onames.append(horigname.split('||')[1])
                else:
                    onames.append(horigname)
        # add ortholg sequence to sequence array
        onames = np.array(onames)
        for e, ename in enumerate(expnames):
            #print ename
            if ename in onames:
                lenorsig = len(sequence[e])
                origseq = np.copy(sequence[e])
                sequence[e] = [sequence[e]]
                for t in np.where(onames == ename)[0]:
                    if len(odomseqs[t]) < lenorsig:
                        print 'too short', len(odomseqs[t]), lenorsig
                        newseqs = np.copy(origseq)
                        off = 0
                        for ll in range(len(odomseqs[t])):
                            if len(odomseqs[t][ll-off]) == len(newseqs[ll]):
                                newseqs[ll] = odomseqs[t][ll-off]
                            else:
                                off += 1
                            #print off
                        odomseqs[t] = np.copy(newseqs)
                        
                    elif len(odomseqs[t]) > lenorsig:
                        print 'too long', len(odomseqs[t]), lenorsig
                        odomseqs[t] = odomseqs[t][:lenorsig]
                    for o, oseqs in enumerate(origseq):
                        if len(oseqs) != len(odomseqs[t][o]):
                            print 'Sequences of orthologs do not match', len(odomseqs[t]), lenorsig, ename, o, len(oseqs), len(odomseqs[t][o]), oseqs, odomseqs[t][o]
                            sys.exit()
                                                 
                    sequence[e].append(odomseqs[t])
                sequence[e] = np.array(sequence[e]).T
            else:
                sequence[e] = np.array([sequence[e]]).T        

elif '--domain_orthologs' in sys.argv:
    domorthfile = sys.argv[sys.argv.index('--domain_orthologs')+1]
    doobj = open(domorthfile, 'r').readlines()
    odnames = []
    odexpnames = []
    odseqs = []
    for d, dline in enumerate(doobj):
        if dline[0] == '>':
            odnames.append(dline[1:].strip())
            horigname = dline[1:].split('||_')[0]
            if '||' in horigname:
                domcoor = dline.split('||_')[0]
                dompi = dposloc(domcoor)
                odexpnames.append([horigname.split('__')[0].split("||")[1], int(domcoor[dompi:])])
            else:
                domcoor = dline.split('||_')[0]
                dompi = dposloc(domcoor)
                odexpnames.append([horigname.split('__')[0], int(domcoor[dompi:])])
            odseqs.append(doobj[d+1].strip())
    
    if len(odnames) == 0:
        print 'No homologs found in file'
        for e, expname in enumerate(expnames):
            for p in range(len(sequence[e])):
                sequence[e][p] = [sequence[e][p]]
    else:
        odexpnames = np.array(odexpnames)
        odseqs = np.array(odseqs)
        if np.sum(np.isin(expnames, odexpnames)) == 0:
            print 'Ortholog names were not compatible'
            sys.exit()
        #print odexpnames
        #print odseqs
        for e, expname in enumerate(expnames):
            potdomseq = np.where(odexpnames[:,0] == expname)[0]
            #print sequence[e]
            #print expname
            #print potdomseq
            if len(potdomseq) == 0:
                for p in range(len(sequence[e])):
                    sequence[e][p] = [sequence[e][p]]
            else:
                potnum = odexpnames[potdomseq,1].astype(int)
                #print odexpnames[potdomseq]
                potdomnum = np.sort(np.unique(potnum))
                #print potdomnum, len(sequence[e])
                #print odseqs[potdomseq]
                #print sequence[e]
                if len(potdomnum) == len(sequence[e]):
                    for p, potdnum in enumerate(potdomnum):
                        if len(sequence[e][p]) == len(odseqs[potdomseq[potnum == potdnum]][0]):
                            sequence[e][p] = np.append([sequence[e][p]],odseqs[potdomseq[potnum == potdnum]])
                        else:
                            print 'Something went wrong, sequences not compatible', expnames[e], p, potdom, len(sequence[e][p]), len(odseqs[potdomseq[potnum == potdnum]][0])
                else:
                    for p, potdnum in enumerate(potdomnum):
                        if len(sequence[e][potdnum]) == len(odseqs[potdomseq[potnum == potdnum]][0]):
                            sequence[e][potdnum] = np.append([sequence[e][potdnum]],odseqs[potdomseq[potnum == potdnum]])
                        else:
                            print 'No homforall: Something went wrong, sequences not compatible', expnames[e], p, potdnum, len(sequence[e][potdnum]),sequence[e][potdnum], len(odseqs[potdomseq[potnum == potdnum]][0]),odseqs[potdomseq[potnum == potdnum]][0]
                            #sys.exit()

            
# remove empty sequences
if '--domainfasta' in sys.argv:
    print 'New sequence cleaning'
    newsequence = []
    for s, seq in enumerate(sequence):
        keep = []
        for t, dseq in enumerate(seq):
            if len(dseq) != 0: 
                keep.append(dseq)
            else:
                print s, t, dseq
        newsequence.append(keep)
    sequence = newsequence


#### reduction of kmer space with reduced alphabet
if '--reduced_alphabet' in sys.argv:
    alphabet = np.genfromtxt(sys.argv[sys.argv.index('--reduced_alphabet')+1], dtype = str)
    aa = np.unique(alphabet[:,1])
    for i, seqarr in enumerate(sequence):
        for j, seq in enumerate(seqarr):
            seq = np.array(list(seq))
            poseq = np.zeros(len(seq), dtype = int)
            for l, oA in enumerate(alphabet[:,0]):
                poseq[seq == oA] = l
            sequence[i][j] = ''.join(alphabet[poseq,1])

else:                            
    aa = ['A', 'R', 'N', 'D', 'P', 'V', 'I', 'C', 'Y', 'H', 'T', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S']

def add_sing(arr, sing):
    outarr = []
    for ar in arr:
        for si in sing:
            outarr.append(ar+si)
    return outarr

##### different possibilities to reduce kmer representation
if '--gappedkmer' in sys.argv:
    gapsize = int(sys.argv[sys.argv.index('--gappedkmer')+1])
    if '--fixcentre' in sys.argv:
        fixc = True
        cent = (k - gapsize)/2
    else:
        fixc = False
    kmers = ['']
    for i in range(k - gapsize):
        kmers = add_sing(kmers, aa)
    kmers = np.array(kmers)
    featmat = np.zeros((len(kmers),len(sequence)), dtype = np.uint8)
    for i, seqarr in enumerate(sequence):
        for seq in seqarr:
            for j in range(len(seq)+1-k):
                fkmer = seq[j:j+k]
                if fixc:
                    featmat[kmers == fkmer[:cent]+fkmer[cent+gapsize:],i] += 1
                else:
                    for l in range(k-gapsize-1):
                        featmat[kmers == fkmer[:1+l]+fkmer[1+l+gapsize:],i] += 1    


elif '--mismatchkmer' in sys.argv:
    gapsize = int(sys.argv[sys.argv.index('--mismatchkmer')+1])
    kmers = ['']
    for i in range(k - gapsize):
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
    allpos = np.arange(k)
    posar = np.arange(gapsize)
    ocheck = True
    while ocheck:
        listpositions.append(np.setdiff1d(allpos,posar))
        posar, ocheck = mix(posar, k-2) # Never get cover the last kmer because otherwise these kmers are counted twice
    listpositions = np.array(listpositions)    
    
    featmat = np.zeros((len(kmers),len(sequence)), dtype = np.uint8)
    for i, seqarr in enumerate(sequence):
        for seq in seqarr:
            for j in range(len(seq)+1-k):
                fkmer = list(seq[j:j+k])
                for lpos in listpositions:
                    mikmer = ''.join(np.array(list(fkmer))[lpos])
                    featmat[kmers == mikmer,i] += 1 

### here only gapsize of possible yet!!!!
elif '--posmismatchkmer' in sys.argv:
    gapsize = int(sys.argv[sys.argv.index('--posmismatchkmer')+1])
    kmerfeat = []
    for kend in range(start-gapsize,k+1-gapsize):    
        kmers = ['']
        for i in range(kend):
            kmers = add_sing(kmers, aa)
        kmerfeat.append(kmers)
    print 'Generate kmer'
    kmerfeat = np.concatenate(kmerfeat)
    rkmerfeat = np.copy(kmerfeat)
    kmerfeat = []
    for g in range(gapsize):
        gappedkmerfeat = []
        for kmer in rkmerfeat:
            gkmer = []
            for i in range(1,len(kmer)):
                gkmer.append(kmer[:i]+'X'+kmer[i:])
            gappedkmerfeat.append(gkmer)
        rkmerfeat = np.unique(np.concatenate(gappedkmerfeat))
    kmerfeat = rkmerfeat
    print 'kmers done ...'

    blsmo = False
    if '--blosumsmooth' in sys.argv:
        blsmo = True
        if '--loadweighttensor' in sys.argv:
            subpos = np.load(sys.argv[sys.argv.index('--loadweighttensor')+1])["subpos"]
            subval = np.load(sys.argv[sys.argv.index('--loadweighttensor')+1])["subval"]
            print kmerfeat[0][0]
            print subpos[0][0]
            print np.array(kmerfeat[0])[subpos[0][0]]
            print subval[0][0]
        else:
            # reads in the blosum matrix file
            blosummatrix = sys.argv[sys.argv.index('--blosumsmooth')+1]
            # defines the maximum number of exchanged amino acids, has to be smaller than  k-1
            maxdiff = int(sys.argv[sys.argv.index('--blosumsmooth')+2])
            
            #background frequencies of each aa in order of array residues
            backfreq = [7.5e-2, 5.2e-2, 4.6e-2, 5.2e-2, 1.8e-2, 4.1e-2, 6.3e-2, 7.1e-2, 2.2e-2, 5.5e-2, 9.1e-2, 5.8e-2, 2.8e-2, 3.9e-2, 5.1e-2, 7.4e-2, 6.0e-2, 1.3e-2, 3.3e-2, 6.5e-2]
            simsmomat = np.genfromtxt(blosummatrix, dtype =str)
            residues = list(simsmomat[:,0])
            simsmomat = simsmomat[:,1:].astype(float)
            
            kmeroutname = 'Residuesmooting'+ str(start)+'-'+str(k) + os.path.splitext(os.path.split(blosummatrix)[-1])[0] + '_maxdiff'+str(maxdiff)+'_gapped'+str(gapsize) 
            # save the prepredictions!!!
            
            subpos, subval, kmeroutname = blosumsub(kmerfeat, simsmomat, maxdiff, kmeroutname)
            #define outname
            np.savez(kmeroutname+'.npz', subpos = subpos, subval = subval)

    if blsmo:
        featmat = np.zeros((len(kmerfeat),len(sequence)), dtype = np.float32)
        def findkmer(seq, kmerset):
            koutarray = np.zeros(len(kmerset))
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
                masks = rmasks

                for p in range(len(seq)+1-l):
                    ckm = np.array(list(seq[p:p+l]))
                    for mask in masks:
                        cckm = np.copy(ckm)
                        cckm[mask] = 'X'
                        ccind = list(kmerset).index(''.join(cckm))
                        koutarray[subpos[ccind]] += subval[ccind]
            return koutarray
    else:
        featmat = np.zeros((len(kmerfeat),len(sequence)), dtype = np.uint8)
        def findkmer(seqar, kmerset, pomask, iy):
            #print seqar
            koutarray = np.zeros(len(kmerset))
            ml = -1
            for l in range(start,k+1):
                ml += 1
                masks = pomask[ml]
                for seq in seqar:     
                    for p in range(len(seq)+1-l):
                        ckm = np.array(list(seq[p:p+l]))
                        #print ckm
                        for mask in masks:
                            cckm = np.copy(ckm)
                            cckm[mask] = 'X'
                            #print cckm
                            koutarray[kmerset == ''.join(cckm)] += 1
            outarrmask = np.where(koutarray > 0)[0]
            koutarray = koutarray[outarrmask]
            #print outarrmask
            #print koutarray
            return outarrmask, koutarray, iy

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

    kmers = np.copy(kmerfeat)
    tmasks = makemask(start,k,gapsize)
    print 'Get k-mers...'
    if mprocessing:
        results = Parallel(n_jobs=num_cores)(delayed(findkmer)(sequence[i], kmers, tmasks, i) for i in range(len(sequence)))
        print 'Distribute multiprocessing results'
        for rpos, res, r in results:
            featmat[rpos, r] = res
    else:
        for i, seqarr in enumerate(sequence):
            rpos, res, r = findkmer(seqarr, kmers, tmasks, i)
            featmat[rpos, r] = res
                  
elif '--seqcomposition' in sys.argv:
    aa = np.sort(aa)
    kmers = []
    for aaa in aa:
        kmers.append(aaa)
    
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

    allpos = np.arange(len(aa))
    lenaa = len(allpos)
    for kend in range(2,k+1):
        posar = np.arange(kend)
        ocheck = True
        while ocheck:
            kmers.append(''.join(aa[posar]))
            posar, ocheck = mix(posar, lenaa-1)
    kmers = np.array(kmers)
    
    featmat = np.zeros((len(kmers),len(sequence)), dtype = np.uint8)
    for i, seqarr in enumerate(sequence):
        for seq in seqarr:
            for j in range(len(seq)+1-k):
                fkmer = list(seq[j:j+k])
                mikmer = ''.join(np.unique(fkmer))
                print fkmer, mikmer
                featmat[kmers == mikmer,i] += 1 

else:    
    kmerfeat = []
    for kend in range(start,k+1):    
        kmers = ['']
        for i in range(kend):
            kmers = add_sing(kmers, aa)
        kmerfeat.append(kmers)
    kmers = np.concatenate(kmerfeat)
# mix position gap with blosum
# mix blosum, homologs and position gap
# introduce N to Blosum: with N == 0 

    if '--Residuesmoothing' in sys.argv:

        if '--loadweighttensor' in sys.argv:
            subpos = np.load(sys.argv[sys.argv.index('--loadweighttensor')+1])["subpos"]
            subval = np.load(sys.argv[sys.argv.index('--loadweighttensor')+1])["subval"]
            print kmerfeat[0][0]
            print subpos[0][0]
            print np.array(kmerfeat[0])[subpos[0][0]]
            print subval[0][0]
        else:
            # reads in the blosum matrix file
            blosummatrix = sys.argv[sys.argv.index('--Residuesmoothing')+1]
            # defines the maximum number of exchanged amino acids, has to be smaller than  k-1
            maxdiff = int(sys.argv[sys.argv.index('--Residuesmoothing')+2])
            
            #background frequencies of each aa in order of array residues
            backfreq = [7.5e-2, 5.2e-2, 4.6e-2, 5.2e-2, 1.8e-2, 4.1e-2, 6.3e-2, 7.1e-2, 2.2e-2, 5.5e-2, 9.1e-2, 5.8e-2, 2.8e-2, 3.9e-2, 5.1e-2, 7.4e-2, 6.0e-2, 1.3e-2, 3.3e-2, 6.5e-2]
            simsmomat = np.genfromtxt(blosummatrix, dtype =str)
            residues = list(simsmomat[:,0])
            simsmomat = simsmomat[:,1:].astype(float)
            
            kmeroutname = 'Residuesmooting'+ str(start)+'-'+str(k) + os.path.splitext(os.path.split(blosummatrix)[-1])[0] + '_maxdiff'+str(maxdiff) 
            subpos, subval, kmeroutname = blosumsub(kmerfeat, simsmomat, maxdiff, kmeroutname)
            #define outname
            # save the prepredictions!!!
            np.savez(kmeroutname+'.npz', subpos = subpos, subval = subval)
        
        print 'Find kmers in protein sequence..'
        ### look for the kmers in the protein sequence and assign weights to them and related kmers. 
        subpos = np.array(subpos)
        subval = np.array(subval)
        kmers = np.concatenate(kmerfeat)
        featmat = np.zeros((len(kmers),len(sequence)))
        for i, seqarr in enumerate(sequence):
            for s, seq in enumerate(seqarr):
                for kend in range(start,k+1):
                    for j in range(len(seq)+1-kend):
                        featmat[np.array(subpos[kmers == seq[j:j+kend]][0]),i] += np.array(subval[kmers == seq[j:j+kend]][0])
                        
    elif '--orthologs' in sys.argv or '--domain_orthologs' in sys.argv:
        if '--orthasone' in sys.argv:
            def orthscore(le):
                return 1
        else:
            def orthscore(le):
                return 1./le
        blosum = False
        pskgap = False
        if '--loadweighttensor' in sys.argv:
            blosum = True
            subpos = np.load(sys.argv[sys.argv.index('--loadweighttensor')+1])["subpos"]
            subval = np.load(sys.argv[sys.argv.index('--loadweighttensor')+1])["subval"]
        
        elif '--posgapkmer' in sys.argv:
            pskgap = True
            gapsize = int(sys.argv[sys.argv.index('--posgapkmer')+1])
            kmerfeat = []
            for kend in range(start-gapsize,k+1-gapsize):    
                kmers = ['']
                for i in range(kend):
                    kmers = add_sing(kmers, aa)
                kmerfeat.append(kmers)
            print 'Generate kmer'
            kmerfeat = np.concatenate(kmerfeat)
            rkmerfeat = np.copy(kmerfeat)
            kmerfeat = []
            for g in range(gapsize):
                gappedkmerfeat = []
                for kmer in rkmerfeat:
                    gkmer = []
                    for i in range(1,len(kmer)):
                        gkmer.append(kmer[:i]+'X'+kmer[i:])
                    gappedkmerfeat.append(gkmer)
                rkmerfeat = np.unique(np.concatenate(gappedkmerfeat))
            mylen = np.vectorize(len)
            rkmerlen = mylen(rkmerfeat)
            rksort = np.argsort(rkmerlen)
            rkmerlen = rkmerlen[rksort]
            kmerfeat = rkmerfeat[rksort]
            for kend in range(start,k+1):
                kmerfeat[rkmerlen == kend] = np.sort(kmerfeat[rkmerlen == kend])
            kmers = np.copy(kmerfeat)
            
            masks = []
            for l in range(start,k+1):
                smasks = []
                for h in range(1,l-1):
                    mask = np.zeros(l) > 0.
                    mask[h] = True
                    smasks.append(mask)
                rmasks = np.copy(smasks)
                for g in range(gapsize-1):
                    smasks = []
                    for mask in rmasks:
                        for h in range(1,l-1):
                            cmask = np.copy(mask)
                            cmask[h] = True
                            if np.sum(cmask) == g + 2:
                                smasks.append(cmask)
                    rmasks = np.unique(smasks, axis = 0) #!!!!
                masks.append(rmasks)

            
            def addkmer(amer, koutarray, kmerset, val, masks):
                ckm = np.array(list(amer))
                for mask in masks:
                    cckm = np.copy(ckm)
                    #print cckm, mask
                    cckm[mask] = 'X'
                    koutarray[kmerset == ''.join(cckm)] += val
        
        elif '--mismatch' in sys.argv:
            pskgap = True
            gapsize = int(sys.argv[sys.argv.index('--mismatch')+1])
            kmerfeat = []
            for kend in range(start-gapsize,k+1-gapsize):    
                kmers = ['']
                for i in range(kend):
                    kmers = add_sing(kmers, aa)
                kmerfeat.append(kmers)
            
            kmerfeat = np.concatenate(kmerfeat)
            kmerfeat = np.array(kmerfeat)
            kmers = np.copy(kmerfeat).astype('S10')
            for f, fkmer in enumerate(kmers):
                kmers[f] = str(int(gapsize))+'$'+fkmer 
            
            masks = []
            for l in range(start,k+1):
                smasks = []
                for h in range(1,l-1):
                    mask = np.ones(l) > 0.
                    mask[h] = False
                    smasks.append(mask)
                rmasks = np.copy(smasks)
                for g in range(gapsize-1):
                    smasks = []
                    for mask in rmasks:
                        for h in range(1,l-1):
                            cmask = np.copy(mask)
                            cmask[h] = False
                            if np.sum(cmask) == l-2-g:
                                smasks.append(cmask)
                    rmasks = np.unique(smasks, axis = 0) 
                masks.append(rmasks)
            

            def addkmer(amer, koutarray, kmerset, val, masks):
                ckm = np.array(list(amer))
                for mask in masks:
                    cckm = np.copy(ckm)
                    #print cckm[mask]
                    koutarray[kmerset == ''.join(cckm[mask])] += val
        
        elif '--gap' in sys.argv:
            pskgap = True
            gapsize = int(sys.argv[sys.argv.index('--gap')+1])
            kmerfeat = []
            for kend in range(start-gapsize,k+1-gapsize):    
                kmers = ['']
                for i in range(kend):
                    kmers = add_sing(kmers, aa)
                kmerfeat.append(kmers)
            
            kmerfeat = np.concatenate(kmerfeat)
            kmers = np.copy(kmerfeat).astype('S10')
            for f, fkmer in enumerate(kmers):
                kmers[f] = str(int(gapsize))+'^'+fkmer        
            
            masks = []
            for l in range(start,k+1):
                smasks = []
                for h in range(1,l-gapsize):
                    mask = np.ones(l) > 0.
                    for hi in range(gapsize):
                        mask[h + hi] = False
                    smasks.append(mask)
                masks.append(smasks)
            #print masks
            
            
            def addkmer(amer, koutarray, kmerset, val, masks):
                ckm = np.array(list(amer))
                for mask in masks:
                    cckm = np.copy(ckm)
                    koutarray[kmerset == ''.join(cckm[mask])] += val
            
        
        
        featmat = np.zeros((len(kmers),len(sequence)))
        # split sequences into lists of amino acids
        print np.shape(featmat)
        from joblib import Parallel, delayed
        import multiprocessing
        if '--mprocessing' in sys.argv:
            num_cores = int(sys.argv[sys.argv.index('--mprocessing')+1])
        else:
            num_cores = 1
        #num_cores = multiprocessing.cpu_count()
        print "Multiprocessing. Using", num_cores, 'cores'
        print 'Find kmers in protein sequence and its homologs..'
        
        def countkmer(seqarr, iy):
            kcountvec = np.zeros(len(kmers))
            #print seqarr
            for s, seq in enumerate(seqarr):
                for kend in range(start,k+1):
                    #print kend
                    domkmer = []
                    for j in range(len(seq[0])+1-kend):
                        jkmer = []
                        for o, orth in enumerate(seq):
                            jkmer.append(orth[j:j+kend])
                        domkmer.append(jkmer)
                    for d in range(len(domkmer)):
                        kat = np.unique(domkmer[d])
                        if len(kat) > 1:
                            delar = []
                            for cs, ks in enumerate(kat):
                                if '?' in ks or '-' in ks or 'X' in ks:
                                    delar.append(cs)
                            kat = np.delete(kat, delar)       
                            #print kat
                        if blosum:
                            for ks in kat:
                                kmind = list(kmers).index(ks)
                                kcountvec[subpos[kmind]] += orthscore(len(kat))*np.array(subval[kmind])
                        elif pskgap:
                            kmasks = masks[kend-start]
                            for ks in kat:
                                addkmer(ks, kcountvec, kmers, orthscore(len(kat)), kmasks)
                        else:
                            for ks in kat:
                                kcountvec[kmers == ks] += orthscore(len(kat))
            return kcountvec, iy
        
        
        results = Parallel(n_jobs=num_cores)(delayed(countkmer)(sequence[i], i) for i in range(len(sequence)))
        for res, r in results:
            featmat[:, r] = res
            
        '''
        for i, seqarr in enumerate(sequence):
            print expnames[i]
            for s, seq in enumerate(seqarr):
                for kend in range(start,k+1):
                    #print kend
                    domkmer = []
                    for j in range(len(seq[0])+1-kend):
                        jkmer = []
                        for o, orth in enumerate(seq):
                            jkmer.append(orth[j:j+kend])
                        domkmer.append(jkmer)
                    for d in range(len(domkmer)):
                        kat = np.unique(domkmer[d])
                        if len(kat) > 1:
                            delar = []
                            for cs, ks in enumerate(kat):
                                if '?' in ks or '-' in ks or 'X' in ks:
                                    delar.append(cs)
                            kat = np.delete(kat, delar)       
                        
                        if blosum:
                            for ks in kat:
                                kmind = list(kmers).index(ks)
                                featmat[subpos[kmind],i] += orthscore(len(kat))*np.array(subval[kmind])
                        elif pskgap:
                            kmasks = masks[kend-start]
                            #print kmasks
                            for ks in kat:
                                addkmer(ks, featmat, kmerfeat, i, orthscore(len(kat)), kmasks)
                        else:
                            for ks in kat:
                                featmat[kmers == ks,i] += orthscore(len(kat))
                '''
    
    else:
        kmers = np.concatenate(kmerfeat)
        #### search for kmers
        if '--sequencescores' in sys.argv:
            featmat = np.zeros((len(kmers),len(sequence)))
        else:
            featmat = np.zeros((len(kmers),len(sequence)), dtype = np.uint8)
        for i, seqarr in enumerate(sequence):
            for s, seq in enumerate(seqarr):
                for kend in range(start,k+1): 
                    for j in range(len(seq)+1-kend):
                        featmat[kmers == seq[j:j+kend],i] += outscore(score[i][s][j:j+kend], scoretype)

        
kmers = kmers[np.sum(featmat, axis = 1)>0]
featmat = featmat[np.sum(featmat, axis = 1)>0]
if '--savetxt' in sys.argv:
    #print featmat
    #for p in range(len(sequence)):
        #print ''.join(sequence[p])
        #for ck in kmers[np.nonzero(featmat[:, p])]:
            #if ck not in ''.join(sequence[p]):
                #print ck
    np.savetxt(outname+"_"+str(k)+"mer_features.txt", np.append([kmers], featmat.T.astype(str), axis = 0).T, fmt = "%s", header = ' '.join(protnames)+"\n"+" ".join(expnames))

elif '--saveorigformat' in sys.argv:
    for i, expn in enumerate(expnames):
        expnames[i] = '"'+expn+'"'
    for i, expn in enumerate(kmers):
        kmers[i] = '"'+expn+'"'
    np.savetxt(outname+"_"+"features.csv", np.append([expnames], featmat.astype(str), axis = 0).T, fmt = "%s", delimiter=",", header = '"",'+','.join(kmers), comments = '')

elif '--savesparse' in sys.argv:
    import scipy.sparse
    print 'Saved in sparse format'
    sparse_matrix = scipy.sparse.csc_matrix(featmat)
    np.savez_compressed(outname+"_"+str(k)+"sparse_mer_features.npz", features=sparse_matrix.data, featshape = sparse_matrix.shape, featindexes = sparse_matrix.indices, featindpr = sparse_matrix.indptr, kmers=kmers, protnames=protnames, expnames=expnames)

else:
    print 'Saved as', outname+"_"+str(k)+"mer_features.npz"
    #print len(np.sum(featmat, axis = 0))
    #print np.where(np.sum(featmat, axis = 0) == 0)[0]
    np.savez_compressed(outname+"_"+str(k)+"mer_features.npz", features=featmat, kmers=kmers, protnames=protnames, expnames=expnames)
