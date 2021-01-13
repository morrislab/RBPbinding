#calculates the pairwise identity between a list of sequences
import numpy as np
import sys, os
from Bio import pairwise2
import operator
import math 
from joblib import Parallel, delayed
import multiprocessing

def readinfasta(seqfile):
    # read in sequences
    gobj = open(seqfile, 'r')
    glines = gobj.readlines()
    i=-1
    starts = []
    totallen = len(glines)
    for gline in glines:
        i+=1
        if ">" in gline:
            starts.append(i)
    #collecting seq names
    seqnames = []
    for sta in starts:
        seqnames.append(glines[sta].split()[0][1:])
    
    sequenz = []
    starts.append(totallen)
    for i in range(len(starts)-1):
        seq = ''
        for j in range(starts[i]+1, starts[i+1]):
            seq = seq + glines[j].strip()
        sequenz.append(seq)
    gobj.close()
    return seqnames, sequenz



if '--mprocessing' in sys.argv:
    num_cores = int(sys.argv[sys.argv.index('--mprocessing')+1])
else:
    num_cores = 1
#num_cores = multiprocessing.cpu_count()


if "--sequences" in sys.argv:
    seqfile = sys.argv[sys.argv.index("--sequences")+1]
    seqnames, sequenz = readinfasta(seqfile)

else:
    print "Please give sequences in form of Fasta file to be aligned"
    sys.exit()



# define the denominator for the sequence identity
# along the alignment, the shorter sequence or the longer sequence
if '--lengthnorm' in sys.argv:
    lnorm = sys.argv[sys.argv.index('--lengthnorm')+1] # min, max, alignment
    def norm(seq1, seq2, al):
        l1 = len(seq1)
        l2 = len(seq2)
        l3 = len(al)
        if lnorm == 'min':
            no = min(l1, l2)
        elif lnorm == 'max':
            no = max(l1,l2)
        elif lnorm == 'alignment':
            no = len(al)
        elif lnorm == 'l1':
            no = l1
        elif lnorm == 'l2':
            no = l2
        else:
            print 'lnorm', lnorm, 'not understood'
            sys.exit()
        return float(no)
else:
    print "please define lengthnorm"
    sys.exit()

# define aligment parameters    
if '--alignment' in sys.argv:
    altype = sys.argv[sys.argv.index('--alignment')+1]
    matrix = sys.argv[sys.argv.index('--alignment')+2]
    gappen = float(sys.argv[sys.argv.index('--alignment')+3])
    gapexten = float(sys.argv[sys.argv.index('--alignment')+4])
    if matrix == 'Blosum62':
        # Define sequence scoring matrix
        from Bio.SubsMat import MatrixInfo as matlist
        matrix = matlist.blosum62

            
    elif matrix == "binary":
        match = float(sys.argv[sys.argv.index('--alignment')+5])
        mismatch = float(sys.argv[sys.argv.index('--alignment')+6])
    else:
        sys.exit()
    
    if altype == 'local' and matrix != "binary":
        def align(seq1, seq2):
            alignment = pairwise2.align.localds(seq1, seq2, matrix, gappen, gapexten)
            return alignment
    elif altype == 'local' and matrix == "binary":
        def align(seq1, seq2):
            alignment = pairwise2.align.localms(seq1, seq2, match, mismatch, gappen, gapexten)
            return alignment
    elif altype == 'global':
        def align(seq1, seq2):
            alignment = pairwise2.align.globalds(seq1, seq2, matrix, gappen, gapexten)
            return alignment        
    elif altype == 'global' and matrix == "binary":
        def align(seq1, seq2):
            alignment = pairwise2.align.globalms(seq1, seq2, match, mismatch, gappen, gapexten)
            return alignment
else:
    print "Define the alignment method"
    sys.exit()
        
if '--outname' in sys.argv:
    outname = sys.argv[sys.argv.index('--outname')+1]
else: 
    outname = os.path.splitext(seqfile)[0]+'_lnorm'+lnorm+'_align'+altype+'_id'



twosets = False
if '--secondset' in sys.argv:
    twosets = True
    secondset = sys.argv[sys.argv.index('--secondset')+1]
    sequencesB = sequenz
    namesB = seqnames
    
    namesA, sequencesA = readinfasta(secondset)
    outname += '.on.'+os.path.splittext(os.path.split(secondset)[1])[0]
    
else:
    namesA = seqnames
    namesB = seqnames
    sequencesA = sequenz    
    sequencesB = sequenz


uniqueseq = False
if '--unique_sequences' in sys.argv:
    uniqueseq = True
    origseqsA = np.copy(sequencesA)
    origseqsB = np.copy(sequencesB)
    sequencesA, returnA= np.unique(sequencesA, return_inverse = True)
    sequencesB, returnB = np.unique(sequencesB, return_inverse = True)
    if len(origseqsA)*len(origseqsB) - len(sequencesA)*len(sequencesB) > 25:
        print "Uniqeseqs considered only"
        print len(origseqsA), len(sequencesA)
        print len(origseqsB), len(sequencesB)
    else:
        sequencesA = origseqsA
        sequencesB = origseqsB
        uniqueseq = False

# Double check that all sequence use upper cases and don't contain unusual letters
if sys.argv[sys.argv.index('--alignment')+2] == 'Blosum62':
    potentiala = np.unique(matrix.keys())
    print 'Clean sequences'
    for s, se in enumerate(sequencesA):
        if se.islower():
            print se
            se = se.upper()
            print se        
        seqarr = np.array(list(se))
        smask = np.isin(seqarr, potentiala)
        if np.sum(~smask) > 0:
            print se
            se = ''.join(seqarr[smask])
            print se
        sequencesA[s] = se
    for s, se in enumerate(sequencesB):
        if se.islower():
            print se
            se = se.upper()
            print se
        seqarr = np.array(list(se))
        smask = np.isin(seqarr, potentiala)
        if np.sum(~smask) > 0:
            print se
            se = ''.join(seqarr[smask])
            print se
        sequencesB[s] = se


# initialize identity matrix
idents = np.ones((len(sequencesB), len(sequencesA)), dtype = np.float16)
print "RUNNING"
if twosets:
    def idtoset(seqA, seqsB, la):
        print la
        idarray = np.zeros(len(seqsB))
        for m in range(len(seqsB)):
            seqB = seqsB[m]
            aligned = align(seqA, seqB)
            idarray[m] = map(operator.eq, aligned[0][0], aligned[0][1]).count(True)/norm(seqA, seqB, aligned[0][0])    
        return la, idarray
    
    if num_cores > 1:
        results = Parallel(n_jobs=num_cores)(delayed(idtoset)(sequencesA[l], sequencesB, l) for l in range(len(sequencesA)))
        for la, isarr in results:
            idents[:, la] = isarr 
    else:
        for l in range(len(sequencesA)):
            l, idarray = idtoset(sequencesA[l], sequencesB, l)
            idents[:, l] = idarray
            
else:
    print 'single set'
    def idtoset(seqA, seqsB, la):
        print la
        idarray = np.zeros(len(seqsB))
        for m in range(la+1, len(seqsB)):
            seqB = seqsB[m]
            aligned = align(seqA, seqB)
            idarray[m] = map(operator.eq, aligned[0][0], aligned[0][1]).count(True)/norm(seqA, seqB, aligned[0][0])    
        return la, idarray
    
    if num_cores > 1:
        results = Parallel(n_jobs=num_cores)(delayed(idtoset)(sequencesA[l], sequencesB, l) for l in range(len(sequencesA)))
        for la, isarr in results:
            idents[:, la] = isarr 
    else:
        for l in range(len(sequencesA)):
            l, idarray = idtoset(sequencesA[l], sequencesB, l)
            idents[:, l] = idarray
    
    
if twosets == False:
    for l in range (len(sequencesA)):
        for m in range(l+1,len(sequencesB)):
            idents[l,m] = idents[m,l]

np.fill_diagonal(idents, 1.)
idents = np.around(idents*100., 2)


if uniqueseq:
    print 'Return to original'
    print np.shape(idents)
    fidents = idents[returnB][:, returnA]
    idents = fidents
    print np.shape(idents)


np.savez_compressed(os.path.splitext(outname)[0]+'.npz', identmat=idents, names=namesA, names2 = namesB)
if '--savetxt' in sys.argv:
    np.savetxt(outname+'.txt', idents, fmt='%3.2f', header = ' '.join(namesA))
if '--savemax' in sys.argv:
    np.savetxt(os.path.splitext(outname)[0]+'-max.txt', np.array([namesB, np.sort(idents, axis = 1)[:,1]], dtype = str).T, fmt = '%s')




