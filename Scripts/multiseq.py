import numpy as np
import sys, os

fastafile = open(sys.argv[1], 'r').readlines()
homologfile = open(sys.argv[2], 'r').readlines()

fname = []
fseq = []
for l, line in enumerate(fastafile):
    if line[0] == '>':
        fname.append(line[1:].strip())
        fseq.append(fastafile[l+1].strip())

hname = []
hseq = []
name = []
for l, line in enumerate(homologfile):
    if line[0] == '>':
        hname.append(line[1:].strip())
        hseq.append(homologfile[l+1].strip())
        name.append(line[1:].split('||_')[0])
name = np.array(name)
hseq = np.array(hseq)
hname = np.array(hname)

for n, na in enumerate(fname):
    wobj = open(na+'_'+os.path.split(sys.argv[2])[1], 'w')
    wobj.write('>'+na+'\n'+fseq[n]+'\n')
    hloc = np.where(name == na)[0]
    for h in hloc:
        wobj.write('>'+hname[h]+'\n'+hseq[h]+'\n')
    wobj.close()

