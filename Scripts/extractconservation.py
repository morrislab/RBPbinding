import numpy as np
import sys, os
import glob


consfile = glob.glob('*'+sys.argv[1])

names = []
conservation = []
sequences = []
for c, cfile in enumerate(consfile):
    obj = open(cfile, 'r').readlines()
    seq = ''
    cons = []
    for l, line in enumerate(obj):
        line = line.strip().split()
        if line[1] == 'reference' and line[2] == 'sequence:':
            names.append(line[-1])
        if line[0].isdigit():
            seq += line[1]
            cons.append(float(line[-1]))
    conservation.append(np.array(cons))
    sequences.append(seq)
#print names
#print sequences
#print conservation

wobj = open('Conservation'+os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'.txt', 'w')
for n, name in enumerate(names):
    conservation[n][conservation[n] == -1000] = np.amin(conservation[n][conservation[n] != -1000])
    wobj.write('>'+name+'\n'+sequences[n]+'\n'+','.join(np.array(conservation[n],dtype = str))+'\n')





