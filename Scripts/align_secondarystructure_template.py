import numpy as np
import sys, os
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62


ssfile = open(sys.argv[1], 'r').readlines() #Rncmpt.RRMKH.ss
secsequences = []
secstr = []
secnames = []
secpronames = []
for s, sline in enumerate(ssfile):
    if sline[0] == '>':
        secpronames.append(sline[1:].strip())
        secnames.append(sline[1:].strip().split('__')[1])
        secsequences.append(ssfile[s+1].strip())
        secstr.append(ssfile[s+2].strip())
secsequences = np.array(secsequences)
secstr = np.array(secstr)
secnames = np.array(secnames)
secpronames = np.array(secpronames)


fastafile = np.genfromtxt(sys.argv[2], dtype = str)
fnames = []
fseq = []
ftnames = []

for l, line in enumerate(fastafile):
    if line[0] == '>':
        fnames.append(line[1:].strip())
        fseq.append(fastafile[l+1].strip())


idfile = np.genfromtxt(sys.argv[3], dtype = str)
for i, idn in enumerate(idfile):
    idfile[i][1] = idn[1].split('||')[0]+'__'+idn[1].rsplit('__')[-1]


def tempstruc(seq1, seq2, temp):
    struc = ''
    start = np.where(np.array(list(seq1)) != '-')[0][0]
    end = np.where(np.array(list(seq1)) != '-')[0][-1]
    co2 = start-1
    stelement = 'C'
    #print 
    #print seq1
    #print seq2
    
    for s in range(start, end+1):
        if seq2[s] != '-':
            co2 += 1
            if co2 == len(temp):
                stelement = 'C'
            else:
                stelement = temp[co2]
        if seq1[s] != '-':
            struc = struc + stelement
    
    #print seq2.replace('-','')
    #print temp
    #print seq1.replace('-','')
    #print struc
    return struc

obj = open(os.path.splitext(sys.argv[2])[0]+'_secondary_templates.fasta',  'w')
fsecond = []
for f, fname in enumerate(fnames):
    i = list(idfile[:,0]).index(fname)
    besttemp = idfile[i,1]
    s = list(secpronames).index(besttemp)
    align = pairwise2.align.globalds(fseq[f], secsequences[s], matrix, -11, -1)
    fsecond.append(tempstruc(align[0][0], align[0][1], secstr[s]))
    obj.write('>'+fname+'\n'+fseq[f]+'\n'+fsecond[-1]+'\n')


