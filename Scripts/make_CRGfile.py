import numpy as np
import sys, os



obj = open(sys.argv[1], 'r').readlines()
wbj = open(os.path.splitext(sys.argv[1])[0]+'.info.txt', 'w')
for l, line in enumerate(obj):
    if l > 0:
        line = line.strip().split()
        i = int(line[12])
        sp = species[i]
        na = line[1].split('__')
        wbj.write(na[0]+';'+sp+';'+na[1]+';'+line[2]+';'+line[4]+';'+line[5]+';'+line[7]+';'+line[8]+';'+line[9]+';'+str(round(float(line[10]),1))+'\n')
    else:
        species = line.strip('#').strip().split()[13:]
        wbj.write('Protname;Species;Domain;CRMG;edist;closestRNCMPT;Age_CRMG;Age_seq;Parent;ParentSID\n')
        print(species)

