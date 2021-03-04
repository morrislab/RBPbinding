import numpy as np
import sys, os


def reducemot(mot):
    for i in range(len(mot)):
        if mot[i] != 'N':
            st = i
            break
    for i in range(len(mot)):
        if mot[-1-i] != 'N':
            en = len(mot) - i
            break
    return mot[st:en]

pfile = open(sys.argv[1], 'r').readlines()
for l, line in enumerate(pfile):
    line = line.strip()
    if len(line) > 0:
        line = line.split()
        if line[0] == 'Motif':
            print line[1], reducemot(pfile[l+1].split()[1])

