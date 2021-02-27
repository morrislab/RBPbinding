import numpy as np
import sys, os


fastafile = sys.argv[1]
homofile = sys.argv[2]

outname  = os.path.splitext(homofile)[0]

fobj = open(fastafile, 'r').readlines()
hobj = open(homofile, 'r').readlines()

protnames = []
sequences = []
for l, line in enumerate(fobj):
    if line[0] == '>':
        protnames.append(line[1:].strip())
        sequences.append([list(fobj[l+1].strip())])
protnames = np.array(protnames)


hprotnames = []
hsequences = []
for l, line in enumerate(hobj):
    #print line
    if line[0] == '>':
        hprotnames.append(line[1:].strip())
        hsequences.append(list(hobj[l+1].strip()))
hprotnames = np.array(hprotnames)
hsequences = np.array(hsequences)

hprot = []
for h, hprotname in enumerate(hprotnames):
    hprot.append(hprotname.split('||_')[0])
hprot = np.array(hprot)

for p, protname in enumerate(protnames):
    si = np.where(hprot == protname)[0]
    for s in si:
        sequences[p].append(hsequences[s])
    

if '--seqweighting' in sys.argv:
    seqweight = sys.argv[sys.argv.index('--seqweighting')+1] # cutoff, linear
    outname += seqweight
    if seqweight == 'cutoff':
        cutoff = float(sys.argv[sys.argv.index('--seqweighting')+2])
        outname+=str(cutoff)
    else:
        cutoff = None
    
    def identity(seA, seB):
        ids = 0.
        tota = 0.
        for p in range(len(seA)):
            if seA[p] != '-' and seB[p] != '-':
                tota += 1.
                if seA[p] == seB[p]:
                    ids += 1.
        return ids/max(1.,tota)
        
    def createweight(seqarr, seqweight, cutoff):
        idmat = np.ones((len(seqarr), len(seqarr)))
        for s, seqa, in enumerate(seqarr):
            for t in range(s +1, len(seqarr)):
                idmat[s,t] = idmat[t,s] = identity(seqa, seqarr[t])
        if seqweight == 'cutoff':
            idmat[idmat >= cutoff] = 1.
            idmat[idmat < cutoff] = 0.
        return idmat
    
    pseqweight = []
    for p, protname in enumerate(protnames):
        weightmat = createweight(sequences[p], seqweight, cutoff)
        pseqweight.append( 1./(1.+np.sum(weightmat, axis = 1)))
else:
    pseqweight = []
    for p, protname in enumerte(protnames):
        pseqweight.append(np.ones(len(sequences[p])))


if '--countgap' in sys.argv:
    aa = np.sort(['A', 'R', 'N', 'D', 'P', 'V', 'I', 'C', 'Y', 'H', 'T', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', '-'])
    outname += '-gaps'
else:
    aa = np.sort(['A', 'R', 'N', 'D', 'P', 'V', 'I', 'C', 'Y', 'H', 'T', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S'])



pssms = []
pseudo = 0.
for p, protname in enumerate(protnames):
    cseqs = np.array(sequences[p])
    psweights = pseqweight[p]
    meff = np.sum(psweights) #+ pseudo/float(len(psweights))
    #print meff
    psm = np.zeros((len(aa), len(cseqs[0])))+meff/float(len(aa))
    for pos in range(len(cseqs[0])):
        for ai, a in enumerate(aa):
            #print pos, a, np.sum(psweights[cseqs[:,pos] == a]), meff
            psm[ai, pos] += (pseudo + np.sum(psweights[cseqs[:,pos] == a]))/meff
    pssms.append((psm/np.sum(psm, axis = 0)).T)        

    
if '--logistic' in sys.argv:
    outname += '-logtrans'
    def sigmoid(x):
        return 1. / (1. + np.exp(-x))
    #sigmoid = np.vectorize(sigmoid)
    for p, psm in enumerate(pssms):
        pssms[p] = sigmoid(psm)
    
if '--infocontent' in sys.argv:
    outname += '-ictrans'
    bfreq = 1./float(len(aa))
    def ic(x):
        ico = x*np.log(x/bfreq)
        ico[ico<0] = 0.
        return ico
    for p, psm in enumerate(pssms):
        pssms[p] = ic(psm)
        

#outobj = open(outname+'-pssm.txt', 'w')
#for p, protname in enumerate(protnames):
    #outobj.write('>'+protname+'\n')
    #for ai, a in enumerate(aa):
        #outobj.write(a+' '+' '.join(np.around(pssms[p][ai],3).astype(str))+'\n')
    #outobj.write('\n\n')#

refseqs = []
for s, seq in enumerate(sequences):
    refseqs.append(''.join(seq[0]))
np.savez_compressed(outname+'-pssm.npz', pssms = np.array(pssms), aminos = aa, protnames = protnames, sequences = refseqs)


if '--conservation' in sys.argv:
    constype = sys.argv[sys.argv.index('--conservation')+1]
    ## add a way to compute conservation from msa!!!
    if constype == 'entropy':
        def consf(pfm):
            logpfm = np.log(pfm/np.mean(pfm, axis = 0))
            logpfm[logpfm < 0] = 0.
            entropy = np.sum(pfm*logpfm , axis = 1)
            return entropy
        
    elif constype == 'variance':
        def consf(pfm):
            varfreq = np.sum((pfm - np.mean(pfm, axis = 0))**2,axis = 1 )
            return varfreq
    elif constype == 'pairsum':
        matrix = Blosum62
        def consf(pfm):
            sp = np.zeros(len(pfm[0]))
            for l in range(len(pfm[0])):
                spl = 0.
                for i in range(len(aa)):
                    for j in range(i +1, len(aa)):
                        spl += pfm[i,l] * pfm[j,l] * matrx[aa[i], aa[j]]/np.sqrt(matrx[aa[i], aa[i]]*matrx[aa[i], aa[i]])
                sp[l] = spl
            return sp
    protconservation = []
    for p, pssm in enumerate(pssms):
        conservation = consf(pssm)
        protconservation.append(conservation)
        
    outobj = open(outname+'-conservation'+constype+'.txt', 'w')
    for p, protname in enumerate(protnames):
        outobj.write('>'+protname+'\n')
        outobj.write(refseqs[p]+'\n')
        outobj.write(','.join(np.around(protconservation[p],3).astype(str))+'\n')
