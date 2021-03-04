import numpy as np
import sys, os 

fused = open(sys.argv[1], 'r').readlines()
names = []
addition = []
seqs = []
for l, line in enumerate(fused):
    if line[0] == '>':
        names.append(line[1:].rsplit('__',1)[0])
        seqs.append(fused[l+1].strip())
        add = line.rsplit('__',1)[1]
        if '-' in add:
            adds = add.split('-')
            add = ''
            for a in adds:
                add+=a.split('_')[0]+'-'
            add = add.strip('-')
        else:
            add = add.split('_')[0]
        addition.append(add)

names = np.array(names)
seqs = np.array(seqs)
addition = np.array(addition)

out = open(os.path.splitext(sys.argv[1])[0].strip('_fused')+'_combined.fasta', 'w')
uname = np.unique(names)
for u, un in enumerate(uname):
    ids = names == un
    out.write('>'+un+'__'+'-'.join(addition[ids])+'\n'+''.join(seqs[ids])+'\n')





