import numpy as np
import sys, os


confidence = np.genfromtxt(sys.argv[1], dtype = str)
confident = confidence[confidence[:,1].astype(float)<=float(sys.argv[3]),0]
print len(confident), sys.argv[3]

pwmfile = open(sys.argv[2], 'r').readlines()

outfile = open(sys.argv[2].rsplit('.',1)[0]+'_confidence'+sys.argv[3]+'.hmot', 'w')
keep = False
for l, line in enumerate(pwmfile):
    if line[:5] == 'Motif':
        if line.strip().split('\t')[1] in confident:
            keep = True
    	else:
            keep = False
    if keep:    
        outfile.write(line)
 




