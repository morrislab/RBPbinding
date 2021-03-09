import numpy as np
import sys, os

#clfile = np.load(sys.argv[1], allow_pickle = True)
#clusters = clfile['clusters']
#names = clfile['names'].astype(str)
def readpwm(fi):
    fi = open(fi,'r').readlines()
    names =[]
    clusternames = []
    for l, line in enumerate(fi):
        if line[0] == '#':
            line = line.strip().split()
            names.append(line[-1].split(','))
            clusternames.append(line[1])
    return names, clusternames

names, clusternames = readpwm(sys.argv[1])


#cls = np.unique(clusters)
cnames = np.genfromtxt(sys.argv[2], dtype = str)
pnames = cnames[:,1]
rnames = []
for c, cname in enumerate(cnames):
    rnames.append(cnames[c,3]+'('+cnames[c,4]+')')
rnames = np.array(rnames)

for c, cl in enumerate(clusternames):
    #rna = names[clusters == cl]
    print '\nCluster',cl, len(names[c])
    print '\n'.join(rnames[np.isin(pnames, names[c])])

