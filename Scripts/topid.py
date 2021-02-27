import numpy as np
import sys, os

f = np.load(sys.argv[1])
idents = f['identmat']
idents = idents.T
names = f['names']
subnames = f['names2']
for n, name in enumerate(names):
    names[n] = name.split('||')[-1]
for n, name in enumerate(subnames):
    subnames[n] = name.split('__')[0]

outname = os.path.splitext(sys.argv[1])[0]
if '--maxexperiment' in sys.argv:
    maxexp = sys.argv[sys.argv.index('--maxexperiment')+1]
    outname += '-maxexp'+maxexp
    namesnum = []
    for n, name in enumerate(names):
        namesnum.append(int(name[-4:]))
    namesnum = np.array(namesnum)
    keep = namesnum <= int(maxexp[-4:])
    names = names[keep]
    idents = idents[keep]

best = np.argmax(idents, axis = 0)


outobj = open(outname+'_maxid.txt', 'w')
for s in range(len(subnames)):
	outobj.write(subnames[s]+'\t'+names[best[s]]+'\t'+str(np.around(idents[best[s],s],1))+'\n')



