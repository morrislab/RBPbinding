import numpy as np
import sys, os

def readfasta(fasta, struc):
    obj = open(fasta, 'r').readlines()
    sequences = []
    names = []
    intscore = []
    for s, sline in enumerate(obj):
        if sline[0] == '>':
            names.append(sline[1:].strip())
            sequences.append(obj[s+1].strip())
            if struc:
                intscore.append(obj[s+2].strip().split(','))
    if struc:
        return np.array(sequences), np.array(names), np.array(intscore)
    else:
        return np.array(sequences), np.array(names)

def domainlocation(fseq, dseq, start):
    for dl in range(start, len(fseq)):
        if dseq == fseq[dl :dl+ len(dseq)]:
            return dl, dl+len(dseq)

    return start, start+len(dseq)

dseqs, dnames = readfasta(sys.argv[1], False)
for d, dname in enumerate(dnames):
    dnames[d] = dname.split('||')[1].split('__')[0]
udnames = np.unique(dnames)
udseqs = []
for ud in udnames:
    udseqs.append('*'.join(dseqs[dnames == ud]))
dnames = udnames
dseqs = udseqs

seqs, names, interface = readfasta(sys.argv[2], True)


for n, name in enumerate(names):
    d = list(dnames).index(name)
    if '*' not in dseqs[d]:
        print '>'+names[n]
        print seqs[n]
        print ','.join(interface[n])
    else:
        dseqint = np.zeros(len(dseqs[d]))
        dseq = dseqs[d].split("*")
        start = 0
        dst = 0
        for ds in dseq:
            st, en = domainlocation(seqs[n], ds, start)
            start = np.copy(en)
            dseqint[dst:dst+len(ds)] = np.around(np.array(interface[n][st:en],dtype = float),3)
            dst += len(ds)+1
        print '>'+names[n]
        print '*'.join(dseq)
        print ','.join(dseqint.astype(str))




