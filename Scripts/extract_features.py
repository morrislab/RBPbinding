# extract_features.py
import numpy as np
import sys, os
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62


aa = np.sort(['A', 'R', 'N', 'D', 'P', 'V', 'I', 'C', 'Y', 'H', 'T', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S'])
fastafile = sys.argv[1]
outname = os.path.splitext(fastafile)[0]
fobj = open(fastafile, 'r').readlines()

protnames = []
sequences = []
for l, line in enumerate(fobj):
    if line[0] == '>':
        protnames.append(line[1:].strip())
        sequences.append(list(fobj[l+1].strip()))
protnames = np.array(protnames)


features = []
for p, protname in enumerate(protnames):
    features.append(np.zeros((len(sequences[p]),1)))
featcollect = []

if '--onehot' in sys.argv:
    outname += '-onehot'
    featcollect.append(aa)
    def makeonehot(seq):
        onho = np.zeros((len(seq), len(aa)))
        for s, sa in enumerate(seq):
            if sa in aa:
                onho[s,list(aa).index(sa)] = 1.
        return onho
            
    for p, protname in enumerate(protnames):
        onehot = makeonehot(sequences[p])
        features[p] = np.append(features[p], onehot, axis = 1)
   
if '--physicochemico' in sys.argv:
    outname += '-physchem'
    # From: "Identification of 14-3-3 Proteins Phosphopeptide-Binding Specificity Using an Affinity-Based Computational Approach"
    phfeatures = ['hydrophob','hydrophil','hydbond','sidevol','polarity','polarize','sasa','netcharge','mass']
    featcollect.append(phfeatures)
    phyprop = np.array([['A',0.62,-0.5,2,27.5,8.1,0.046,1.181,0.007187,71.0788],
                        ['C',0.29,-1,2,44.6,5.5,0.128,1.461,-0.03661, 103.1388],
                        ['D',-0.9,3,4,40,13,0.105,1.587,-0.02382,115.0886], 
                        ['E',-0.74,3,4,62,12.3,0.151,1.862,0.006802,129.1155],
                        ['F',1.19,-2.5,2,115.5,5.2,0.29,2.228,0.037552,147.1766],
                        ['G',0.48,0,2,0,9,0,0.881,0.179052,57.0519],
                        ['H',-0.4,-0.5,4,79,10.4,0.23,2.025,-0.01069,137.1411],
                        ['I',1.38,-1.8,2,93.5,5.2,0.186,1.81,0.021631,113.1594],
                        ['K',-1.5,3,2,100,11.3,0.219,2.258,0.017708,128.1741],
                        ['L',1.06,-1.8,2,93.5,4.9,0.186,1.931,0.051672,113.1594],
                        ['M',0.64,-1.3,2,94.1,5.7,0.221,2.034,0.002683,131.1986],
                        ['N',-0.78,2,4,58.7,11.6,0.134,1.655,0.005392,114.1039],
                        ['P',0.12,0,2,41.9,8,0.131,1.468,0.239531,97.1167],
                        ['Q',-0.85,0.2,4,80.7,10.5,0.18,1.932,0.049211,128.1307],
                        ['R',-2.53,3,4,105,10.5,0.18,1.932,0.049211,156.1875],
                        ['S',-0.18,0.3,4,29.3,9.2,0.062,1.298,0.004627,87.0782],
                        ['T',-0.05,-0.4,4,51.3,8.6,0.108,1.525,0.003352,101.1051],
                        ['V',1.08,-1.5,2,71.5,5.9,0.14,1.645,0.057004,99.1326],
                        ['W',0.81,-3.4,3,145.5,5.4,0.409,2.663,0.037977,186.2132],
                        ['Y',0.26,-2.3,3,117.3,6.2,0.298,2.368,0.023599,163.1760]]).T
    def makephysico(seq):
        onho = np.zeros((len(seq), len(phfeatures)))
        for s, sa in enumerate(seq):
            if sa in phyprop[0]:
                onho[s,:] = phyprop[1:, phyprop[0] == sa].T[0].astype(float)
        return onho
    
    for p, protname in enumerate(protnames):
        phys = makephysico(sequences[p])
        features[p] = np.append(features[p], phys, axis = 1)


if '--pssm' in sys.argv:
    pssmfile = sys.argv[sys.argv.index('--pssm')+1]
    outname += '-PSSM'+os.path.splitext(os.path.split(pssmfile)[1])[0]
    pfile = np.load(pssmfile, allow_pickle = True)
    psname = pfile['protnames']
    pssms = pfile['pssms']
    aminos = pfile['aminos']
    pamino = []
    for ai, ab in enumerate(aminos):
        pamino.append('pssm'+ab)
    for p, protname in enumerate(protnames):
        pf = list(psname).index(protname)
        features[p] = np.append(features[p], pssms[pf], axis = 1)
    featcollect.append(pamino)

if '--disordered' in sys.argv:
    outname += '-disord'
    disfile = sys.argv[sys.argv.index('--disordered')+1]
    pfile = open(disfile, 'r').readlines()
    disname = []
    disorder = []
    fesequence = []
    for l, line in enumerate(pfile):
        if line[0] == '>':
            disname.append(line[1:].strip())
            fesequence.append(pfile[l+1].strip())
        if line.strip() == '#diso':
            disord = np.array(pfile[l+1].strip().split(','), dtype = float)
            disorder.append(disord.reshape(len(disord),1))
    
    for p, protname in enumerate(protnames):
        #print protname
        if protname in disname:
            pf = list(disname).index(protname)
            if len(sequences[p]) != fesequence[pf]:
                #print sequences[p], len(fesequence[pf])
                alignment = pairwise2.align.globalds(fesequence[pf], ''.join(sequences[p]), matrix, -11, -1)
                #print alignment[0][0]
                #print alignment[0][1]
                alignmask = np.array(list(alignment[0][0])) != '-'
                newdis = np.zeros(len(sequences[p]))
                newdis[alignmask] = disorder[pf].T[0]
                newdis = np.array([newdis]).T
                disorder[pf] = newdis
        else:
            newdis = np.zeros((len(sequences[p]),1))
        features[p] = np.append(features[p], newdis, axis = 1)

    featcollect.append(['DISO'])

if '--acc' in sys.argv:
    outname += '-acc'
    asafile = sys.argv[sys.argv.index('--acc')+1]
    pfile = open(asafile, 'r').readlines()
    asaname = []
    asas = []
    fesequence = []
    for l, line in enumerate(pfile):
        if line[0] == '>':
            asaname.append(line[1:].strip())
            fesequence.append(pfile[l+1].strip())
        if line.strip() == '#acc':
            asa = []
            for i in range(3):
                asa.append(pfile[l+1+i].strip().split(','))
            asas.append(np.array(asa).T)
    
    for p, protname in enumerate(protnames):
        if protname in asaname:
            pf = list(asaname).index(protname)
            if len(sequences[p]) != fesequence[pf]:
                alignment = pairwise2.align.globalds(fesequence[pf], ''.join(sequences[p]), matrix, -11, -1)
                alignmask = np.array(list(alignment[0][0])) != '-'
                newdis = np.zeros((len(sequences[p]),3))
                newdis[alignmask] = asas[pf]
                asas[pf] = newdis
        else:
            newdis = np.zeros((len(sequences[p]),3))
        features[p] = np.append(features[p], newdis, axis = 1)
    asofeat = []
    for i in range(3):
        asofeat.append('ACC'+str(i))
    featcollect.append(asofeat)

    
if '--secstructure' in sys.argv:
    outname += '-2ndstr'
    secfile = sys.argv[sys.argv.index('--secstructure')+1]
    sectype = sys.argv[sys.argv.index('--secstructure')+2]
    outname += sectype
    pfile = open(secfile, 'r').readlines()
    secname = []
    secs = []
    fesequence = []
    readsec = False
    for l, line in enumerate(pfile):
        if line[0] == '>':
            secname.append(line[1:].strip())
            fesequence.append(pfile[l+1].strip())
        if line.strip() == '#'+sectype:
            readsec = True
            sec = []
        elif readsec and line[0] != '#':
            sec.append(line.strip().split(','))
        elif readsec and (line[0] == '#' or line[0] == '>'):
            #print sec
            secs.append(np.array(sec, dtype = float).T)
            readsec = False
            
        
    for p, protname in enumerate(protnames):
        if protname in secname:
            pf = list(secname).index(protname)
            if len(sequences[p]) != fesequence[pf]:
                alignment = pairwise2.align.globalds(fesequence[pf], ''.join(sequences[p]), matrix, -11, -1)
                alignmask = np.array(list(alignment[0][0])) != '-'
                newdis = np.zeros((len(sequences[p]),len(secs[0][0])))
                newdis[alignmask] = secs[pf]
                secs[pf] = newdis
        else:
            newdis = np.zeros((len(sequences[p]),len(secs[0][0])))
        features[p] = np.append(features[p], newdis, axis = 1)
    
    ssfeat = []
    for i in range(len(secs[0][0])):
        ssfeat.append('SS'+str(i))
    featcollect.append(ssfeat)

    
if '--conservation' in sys.argv:
    outname += '-conservation'
    consfile = sys.argv[sys.argv.index('--conservation')+1]
    pfile = open(consfile, 'r').readlines()
    consname = []
    conservation = []
    fesequence = []
    for l, line in enumerate(pfile):
        if line[0] == '>':
            consname.append(line[1:].strip())
            fesequence.append(pfile[l+1].strip())
            conser = np.array(pfile[l+2].strip().split(','), dtype = float)
            conservation.append(conser.reshape(len(conser),1))
    
    for p, protname in enumerate(protnames):
        if protname in consname:
            pf = list(consname).index(protname)
            if len(sequences[p]) != fesequence[pf]:
                alignment = pairwise2.align.globalds(fesequence[pf], ''.join(sequences[p]), matrix, -11, -1)
                alignmask = np.array(list(alignment[0][0])) != '-'
                newdis = np.zeros((len(sequences[p]),len(conservation[0][0])))
                newdis[alignmask] = conservation[pf]
                conservation[pf] = newdis
        else:
            newdis = np.zeros((len(sequences[p]),len(conservation[0][0])))
        features[p] = np.append(features[p], newdis, axis = 1)
    featcollect.append(['Conser'])


for f, feat in enumerate(features):
    features[f] = feat[:,1:]
featcollect = np.concatenate(featcollect)
print featcollect
outsequences = []
for seq in sequences:
    outsequences.append(''.join(seq))


np.savez_compressed(outname+'-features.npz', featmat = features, features = featcollect, protnames = protnames, sequences = outsequences)
    
    
