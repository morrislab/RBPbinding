import numpy as np
import sys, os
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

scorefasta = open(sys.argv[1], 'r').readlines()
scoreline = int(sys.argv[2])

fnames = []
fscore = []
fseqs = []

for l, line in enumerate(scorefasta):
    if line[0] == '>':
        fnames.append(line[1:].strip())
        fseq = scorefasta[l+1].strip()
        if '*' in fseq:
            fseq = np.array(list(fseq))
            fmask = fseq != '*'
            fseq = ''.join(fseq[fmask])
        else:
            fmask = np.ones(len(fseq)) == 1
        fscore.append(np.array(scorefasta[l+scoreline].strip().split(','))[fmask])
        fseqs.append(fseq)
fnames = np.array(fnames)
fscore = np.array(fscore)
fseqs = np.array(fseqs)
if '--domainfasta' in sys.argv:
    fsort = np.argsort(fnames)
    fnames = fnames[fsort]
    fscore = fscore[fsort]
    fseqs = fseqs[fsort]
    for f, fname in enumerate(fnames):
        fnames[f] = fname.split('__')[0]
    ufnames = np.unique(fnames)
    ufseqs = []
    ufscore = []
    for f, fname in enumerate(ufnames):
        ufseqs.append(''.join(fseqs[fnames == fname]))
        ufscore.append(np.concatenate(fscore[fnames == fname]))
    fnames = ufnames
    fscore = np.array(ufscore)
    fseqs = np.array(ufseqs)


intfasta = open(sys.argv[3], 'r').readlines()

inames = []
iseqs = []
interface = []
idistance = []
for l, line in enumerate(intfasta):
    if line[0] == '>':
        inames.append(line[1:].strip().split(':'))
        iseqs.append(intfasta[l+1].strip())
        interface.append(intfasta[l+2].strip())




outobj = open(os.path.splitext(sys.argv[1])[0]+os.path.splitext(os.path.split(sys.argv[3])[1])[0]+'.txt','w')
for i, iname in enumerate(inames):
    found = False
    for icont in iname:
        if icont in fnames:
            f = list(fnames).index(icont)
            found = True
            break
    if found:
        imask = np.array(list(interface[i]))!='X'
        iseq = ''.join(np.array(list(iseqs[i]))[imask])
        fseq = fseqs[f]
        alignment = pairwise2.align.globalds(iseq, fseq, matrix, -22,-2)
        k = 0
        iscore = []
        for j, ie in enumerate(alignment[0][0]):
            if ie !='-':
                if alignment[0][1][j] != '-':
                    iscore.append(fscore[f][k])
                    k+=1
                else:
                    iscore.append('?')
            elif alignment[0][1][j] != '-':
                k +=1
        #print iname, fnames[f]
        #print ''.join(np.array(list(alignment[0][0]))[np.array(list(alignment[0][0]))!= '-'])
        #print ''.join(np.array(list(alignment[0][1]))[np.array(list(alignment[0][0]))!= '-'])
        #print iscore
        iout = np.array(list(interface[i]), dtype = '|S20')
        iout[imask] = np.array(iscore)
        outobj.write('>'+iname[0]+'\n'+iseqs[i]+'\n'+','.join(iout)+'\n')
    else:
        print '\n',iname, 'not in score file','\n'
        sys.exit()
