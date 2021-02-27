import numpy as np
import sys,os

infofile = open(sys.argv[1], 'r').readlines()

names =[]
domains = []
locations = []
sequences = []
dnums = []

fasts = open(os.path.splitext(sys.argv[1])[0]+'-concatenated.fasta', 'w')

for l, line in enumerate(infofile):
    line = line.strip().split()
    names.append('>'+line[0].strip())
    dnum  = int(line[1])
    dnums.append(dnum)
    domains.append(line[2: 2+dnum])
    locations.append(line[2+dnum:2+dnum + dnum*2])
    sequences.append(line[-dnum:])

for d, doms in enumerate(domains):
    for e, dom in enumerate(doms):
        if not dom[-1].isdigit():
            domains[d][e] = dom+'_'+str(e)

for n, name in enumerate(names):
    if dnums[n] == 1:
        fasts.write(name+'__'+domains[n][0]+'\n')
        fasts.write(sequences[n][0]+'\n')
    elif dnums[n] > 1:
        d = 0
        while True:
            print d, sequences[n], locations[n], domains[n]
            if locations[n][d*2+1] == locations[n][d*2+2]:
                sequences[n][d] = sequences[n][d]+sequences[n][d+1]
                sequences[n].remove(sequences[n][d+1])
                locations[n] = locations[n][:d*2+1]+locations[n][d*2+3:]
                domains[n][d] = domains[n][d]+'-'+domains[n][d+1]
                del domains[n][d+1]
            else:
                d += 1
            if d == len(domains[n]) -1:
                break
            print sequences[n], locations[n], domains[n]
        for d, dom in enumerate(domains[n]):
            fasts.write(name+'__'+dom+'\n')
            fasts.write(sequences[n][d]+'\n')




