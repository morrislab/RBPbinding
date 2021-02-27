import numpy as np
import glob
import sys, os

ending = sys.argv[1]
outname = sys.argv[2]

files = np.sort(glob.glob('*'+ending))

mobj = open(outname + '_motifs'+ending, 'w')
iobj = open(outname +ending, 'w')
for f, fi in enumerate(files):
    obj = open(fi, 'r').readlines()
    print fi
    motif = ''
    name = ''
    for l, line in enumerate(obj):
        if line[0] == '#':
            motif += line[1:].strip()+','
        elif line[0] == '>':
            name = line.strip()
            iobj.write(line)
        else:
            iobj.write(line)
    if motif != '':
        mobj.write(name+'\t'+motif.strip(',')+'\n')



