import numpy as np
import sys, os

tsets = sys.argv[1]
if ',' in tsets:
    tsets = tsets.strip().split(',')
elif '-' in tsets:
    tsets = tsets.strip().split('-')
    tsets = list(np.arange(int(tsets[0]), int(tsets[1])).astype(str))

tsetfile = sys.argv[2]
idfile = sys.argv[3]

if os.path.splitext(idfile)[-1] == '.npz':
    npidfile = np.load(idfile)
    idmat = npidfile['identmat']
    idnames = list(npidfile['names'])
    #print idnames
else:
    idnames = open(idfile, 'r').readline().strip().split()[1:]
    idmat = np.genfromtxt(idfile)
expnames = []
if '||' in idnames[0]:
    for i, iname in enumerate(idnames):
        expnames.append(iname.split("||")[-1])
else:
    expnames = idnames
if '>' in expnames[0]:
    for e, expn in enumerate(expnames):
        expnames[e] = expn[1:]

simil = False
if '--similaritymat' in sys.argv:
    simil = True
    simnames = open(sys.argv[sys.argv.index('--similaritymat')+1], 'r').readline().strip().split()[1:]
    similname = os.path.splitext(os.path.split(sys.argv[sys.argv.index('--similaritymat')+1])[1])[0]
    simat = np.genfromtxt(sys.argv[sys.argv.index('--similaritymat')+1])
    if '||' in simnames[0]:
        for i, iname in enumerate(simnames):
            simnames[i] = iname.split("||")[-1]

    if '>' in simnames[0]:
        for e, expn in enumerate(simnames):
            simnames[e] = expn[1:]





tobj = open(tsetfile, 'r').readlines()
setnum = []
trainsets = []
testsets = []
for l, line in enumerate(tobj):
    if line[:6] == '###Set':
        if line.strip().split()[1] in tsets:
            setnum.append(line.strip().split()[1])
            if tobj[l + 2].split()[0] in expnames:
                train = tobj[l+2].strip().split()
                if tobj[l+3].strip() == '##Test:':
                    test = tobj[l+4].strip().split()
                else:
                    test = tobj[l+5].strip().split()
            else:
                train = tobj[l+3].strip().split()
                test = tobj[l+6].strip().split()
            trainsets.append(train)
            testsets.append(test)

    
if '--topn' in sys.argv:
    topn = int(sys.argv[sys.argv.index('--topn')+1])
else:
    topn = 1

savet = False
if '--savetxt' in sys.argv:
    if len(sys.argv) > sys.argv.index('--savetxt')+1:
        outname = sys.argv[sys.argv.index('--savetxt')+1]
    else:
        outname = ''
    savet = True
    if simil:
        out = open(outname+os.path.splitext(os.path.split(sys.argv[2])[1])[0]+'_'+os.path.splitext(os.path.split(sys.argv[3])[1])[0]+'_bid-'+similname+'.dat', 'w')
    outid = open(outname+os.path.splitext(os.path.split(sys.argv[2])[1])[0]+'_'+os.path.splitext(os.path.split(sys.argv[3])[1])[0]+'_bid.dat', 'w')


if simil:
    # align two sets to each other
    comset = np.intersect1d(expnames, simnames)
    idmask = np.isin(expnames, comset)
    simmask = np.isin(simnames, comset)
    idmat = idmat[idmask][:, idmask]
    simat = simat[simmask][:, simmask]
    expnames = list(np.array(expnames)[idmask])
    simnames = list(np.array(simnames)[simmask])
    
    
for tx, tset in enumerate(testsets):
    train = trainsets[tx]
    trainmat = idmat[:, np.isin(expnames, train)]
    trainexp = np.array(expnames)[np.isin(expnames, train)]
    #print np.shape(trainmat)
    for tpr in tset:
        if tpr in expnames:
            besttrains = np.argsort(trainmat[expnames.index(tpr)])[-topn:]
            bestp = trainexp[besttrains]
            bestrainsid = trainmat[expnames.index(tpr)][besttrains]
            if simil:
                indexsim = []
                for bp in bestp:
                    indexsim.append(simnames.index(bp))
                if savet:
                    out.write(tpr+' '+' '.join(bestp)+' '+' '.join(simat[simnames.index(tpr), indexsim].astype(str))+'\n')
                print tpr+' '+' '.join(bestp)+' '+' '.join(simat[simnames.index(tpr), indexsim].astype(str))
            else:
                print tpr+' '+' '.join(bestp)+' '+' '.join(bestrainsid.astype(str))
            if savet:
                outid.write(tpr+' '+' '.join(bestp)+' '+' '.join(bestrainsid.astype(str))+'\n')
                
        #else:
            #print tpr
        
        
        
    

