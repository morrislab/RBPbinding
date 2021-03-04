import numpy as np
import sys,os 

def keepfunc(names, matrix, rightnames):
    keep = []
    for rightname in rightnames:
        keep.append(list(names).index(rightname))
    names = np.array(names)[keep]
    matrix = matrix[keep]
    matrix = matrix[:,keep]
    return names, matrix

def readin(plotfeatures, featureaxis, prots):
    if os.path.splitext(plotfeatures)[-1] == '.npz':
        obj = np.load(plotfeatures)
        ofiles = obj.files
        if 'features' in ofiles:
            plotfeatures = obj['features']
            pnames = obj['expnames']
        elif 'profiles' in ofiles:
            plotfeatures = obj['profiles']
            pnames = obj['names']
    else:
        pnames = open(plotfeatures, 'r').readline().strip().split()[1:]
        if '||' in pnames[0]:
            for c, cname in enumerate(pnames):
                pnames[c] = cname.split('||')[-1]
        plotfeatures = np.genfromtxt(plotfeatures)
    
    if featureaxis == '0':
        plotfeatures = plotfeatures.T
    
    if prots is not None:
        if featureaxis == 'both':
            pnames, plotfeatures = keepfunc(pnames, plotfeatures, prots)
            np.fill_diagonal(plotfeatures, np.amax(plotfeatures))
        else:
            keep = []
            for cname in prots:
                keep.append(pnames.index(cname))
            pnames = np.array(pnames)[keep]
            plotfeatures = plotfeatures[keep]
                
    prots = pnames[:]

    return pnames, plotfeatures, prots

def check(val):
    if val.upper() == 'TRUE':
        return True
    else: 
        return False

# complete list of proteins that will be split
pset = np.genfromtxt(sys.argv[1], dtype = str)
if len(np.shape(pset)) == 2:
    pset = pset.T
else:
    pset = np.array([pset, pset])
outname = os.path.splitext(sys.argv[1])[0]

### similarity file allows to have only proteins in trainset that have less than a certain similarity
simfile= sys.argv[2]
if simfile != "None":
    simcuts = float(sys.argv[3])
    outname += '-'+os.path.splitext(os.path.split(sys.argv[2])[1])[0]+str(simcuts)
    snames, smat, prots = readin(sys.argv[2], 'both', pset[1])

# idfile contains identities and can be used to filter for proteins with a certain identity in training set
idfile = sys.argv[4]
if idfile != 'None':
    idcuts = float(sys.argv[5])
    outname += '-'+os.path.splitext(sys.argv[4])[0]+str(idcuts)
    idnames, idmat, prots = readin(sys.argv[4], 'both', pset[1])

if simfile != 'None' and idfile != 'None':
    if not np.array_equal(snames, idnames):
        ssort = np.argsort(snames)
        idsort = np.argsort(idnames)
        ssort = ssort[np.isin(snames, idnames)]
        idsort = idsort[np.isin(idnames, snames)]
        snames = snames[ssort]
        idnames = idnames[idnames]
        smat = smat[ssort][:, ssort]
        idmat = idmat[idsort][:, idsort]
        prots = idnames

# if 'True' the proteins will be filtered so that only one protein exists that is close to another protein in sequence identity or specificity
# To avoid high correlation due to lots of similar proteins
inbetween = check(sys.argv[6])

if inbetween: 
    outname += '-inbtwfilter'
    inbmat = np.ones((len(prots), len(prots))) == 1
    if simfile != 'None':
        inbmat *= smat > simcuts
    if idfile != "None":
        inbmat *= idmat > idcuts
    rem = []
    for p in range(len(prots)):
        if p not in rem:
            potrem = np.where(inbmat[p, p+1:])[0]
            if len(potrem) > 0:
                rem = np.append(rem,potrem + p+1)
    rem = np.unique(rem)
    mask = np.delete(np.arange(len(prots)), rem)
    prots = prots[mask]
    if simfile != "None":
        smat = smat[mask][:, mask]
    if idfile != 'None':
        idmat = idmat[mask][:, mask]

if idfile == 'None' and simfile == 'None':
    prots = pset[1]
    prots2 = pset[0]
else:
    prots2 = pset[0][np.isin(pset[1],prots)]



numtest = sys.argv[7]
outname += '_testset'+numtest
if numtest == 'Single':
    numtest = len(prots)
else:
    numtest = int(numtest)

# permut proteins
permut = np.random.permutation(len(prots))
prots = prots[permut]
prots2 = prots2[permut]
if simfile != 'None':
    smat = smat[permut][:, permut]
if idfile != "None":
    idmat = idmat[permut][:, permut]

wobj = open(outname+'.list' , 'w')

testsize = len(prots)/numtest
left = len(prots)%numtest
le = 0


if simfile == "None" and idfile == 'None':
    inbmat = np.zeros((len(prots), len(prots))) == 1
else:
    inbmat = np.ones((len(prots), len(prots))) == 1
if simfile != 'None':
    inbmat *= smat > simcuts
if idfile != "None":
    inbmat *= idmat > idcuts

for p in range(numtest):
    testset = np.arange(p*testsize+le, (p+1)*testsize+le+int(le < left))
    simset = np.where(np.sum(inbmat[testset], axis = 0) > 0)[0]
    if le < left:
        le += 1
    print len(testset), len(simset)
    trainset = np.delete(np.arange(len(prots)), np.unique(np.append(testset, simset)))
    print len(trainset)
    wobj.write('###Set '+str(p)+'\n'+'##Train:'+'\n'+' '.join(prots[trainset])+'\n'+' '.join(prots2[trainset])+'\n'+'##Test:'+'\n'+' '.join(prots[testset])+'\n'+' '.join(prots2[testset])+'\n')









