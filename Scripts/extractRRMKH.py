
import numpy as np
import sys, os


listr=np.genfromtxt(sys.argv[1], dtype = str)

fasta = open(sys.argv[2], 'r').readlines()
for l, line in enumerate(fasta):
    if line[0] == '>':
        if line.split('||')[1].strip().split('_')[0] in listr:
            print line.strip()
            print fasta[l+1].strip()




