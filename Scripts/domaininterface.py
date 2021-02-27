import numpy as np
import sys, os

domainfasta = open(sys.argv[1], 'r').readlines()

fnames = []
fprot = []
fseqs = []

for l, line in enumerate(domainfasta):
    if line[0] == '>':
        fnames.append(line.strip())
        fprot.append(line.strip().split('__')[0])
        fseqs.append(domainfasta[l+1].strip())


intfasta = open(sys.argv[2], 'r').readlines()

inames = []
iseqs = []
iinfo = []
for l, line in enumerate(intfasta):
    if line[0] == '>':
        inames.append(line.strip())
        iseqs.append(''.join(intfasta[l+1].strip().split(',')))
        info = []
        for j in range(6):
            info.append(intfasta[l+1+j].strip().split(','))
        iinfo.append(np.array(info))

def off(lseq, shseq):
    found = False
    for d in range(len(lseq) - len(shseq)+1):
        if lseq[d:len(shseq)+d] == shseq:
            found = True
            break
    if found == False:
        print 'not found'
        print lseq
        print shseq
        sys.exit()
    return d


outobj = open(os.path.splitext(sys.argv[1])[0]+'.int', 'w')
for f, fname in enumerate(fnames):
    outobj.write(fname+'\n')
    i = inames.index(fprot[f])
    offis = off(iseqs[i], fseqs[f])
    endis = offis + len(fseqs[f])
    for j in range(6):
        outobj.write(','.join(iinfo[i][j][offis:endis].astype(str))+'\n')






