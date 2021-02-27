import numpy as np
import sys, os


confidence = np.genfromtxt(sys.argv[1], dtype = str)
confident = confidence[confidence[:,1].astype(float)<=float(sys.argv[3]),0]
print len(confident), sys.argv[3]

zfile = np.genfromtxt(sys.argv[2], dtype = str)
znames = open(sys.argv[2]).readline().strip().split()[1:]
keep = []
for z,zn in enumerate(znames):
    if zn.strip('.Z') in confident:
        keep.append(z)
#print keep, len(keep)
np.savetxt(sys.argv[2].rsplit('.',1)[0]+'_confidence'+sys.argv[3]+'.txt', zfile[:,keep], header = ' '.join(np.array(znames)[keep]), fmt='%s')




