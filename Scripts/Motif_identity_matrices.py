# computes similarities between RNAcompete measurements
# run:
# python ~/Work/RBPs/Scripts/motif_identity/Motif_identity_matrices.py Zscores_426_origrncmpt.txt --zscores --outname Zscores_426_origrncmpt_pearson.mat --pearson --savetxt
#python ~/Work/RBPs/Scripts/motif_identity/Motif_identity_matrices.py Zscores_426_origrncmpt.txt --zscores --outname Zscores_426_origrncmpt_top10.mat --topnum 10 --topset --savetxt
import numpy as np
import sys, os
import glob
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import cosine
from sklearn.metrics import r2_score

zfile = sys.argv[1]


if "--zscores" in sys.argv:
    motnames = open(zfile, 'r').readline().strip().split()[1:]
    #mer7 = np.genfromtxt(zfile, dtype = str)
    #zscores = mer7[:,1:].astype(float)
    #mer7 = mer7[:,0]
    obj = open(zfile, 'r').readlines()[1:]
    mer7 = []
    zscores = []
    for l, line in enumerate(obj):
        line = line.strip().split()
        mer7.append(line[0])
        zscores.append(line[1:])
    mer7 = np.array(mer7)
    zscores = np.array(zscores, dtype = float)

elif "--affprofiles" in sys.argv:
    zscores = np.genfromtxt(zfile, skip_header = 1)[:,3:]
    mer7 = np.genfromtxt(zfile, skip_header = 1, dtype = str)[:,2]
    motnames = np.genfromtxt(zfile, dtype = str)[0, 3:]
elif "--files" in sys.argv:
    filenames = sys.argv[sys.argv.index('--files')+1]
    listfiles = glob.glob('*'+filenames)
    zscores = []
    motnames = []
    for l, lfile in enumerate(listfiles):
        zscores.append(np.genfromtxt(lfile)[:,-1])
        motnames.append(lfile.split('_')[1])
    zscores = np.array(zscores).T
else:
    print "Please provide data type"
    sys.exit()

zscores = np.nan_to_num(zscores)
lenvec = len(zscores)
motnames = np.array(motnames)



if '--negZ' in sys.argv:
    zscores = -zscores

if '--normZ2' in sys.argv:
    print "Z2 norm"
    zscores = zscores/np.sqrt(np.sum(zscores*zscores, axis = 0))


if '--normZ1' in sys.argv:
    zscores = zscores/np.sum(np.absolute(zscores), axis = 0)

if '--preprocessingZ' in sys.argv:
   zscores = (preprocessing.scale(zscores.T)).T


if '--binaryZ' in sys.argv:
    print 'binaryZ'
    #print np.shape(zscores)
    bcoff = float(sys.argv[sys.argv.index('--binaryZ')+1])
    bctype = sys.argv[sys.argv.index('--binaryZ')+2]
    if bctype == 'cut':
        if bcoff < 1:
            bcoff = int(np.shape(zscores)[-1] * bcoff)
        else:
            bcoff = int(bcoff)
        for py, Yp in enumerate(zscores):
            sort = np.argsort(-Yp)
            zscores[py, sort[:bcoff]] = 1
            zscores[py, sort[bcoff:]] = 0
    if bctype == 'sig':
        bcoffset = -np.sort(-zscores.flatten())
        bcoff = bcoffset[int(len(bcoffset)*bcoff)]
        print 'Y lager than', bcoff, np.shape(zscores)
        for py, Yp in enumerate(zscores):
            ymask = Yp>bcoff
            zscores[py, ymask] = 1
            zscores[py, ~ymask] = 0

zscores = np.nan_to_num(zscores)




if "--SVD" in sys.argv:
    perP = float(sys.argv[sys.argv.index("--SVD")+1])
    
    print "SVD\nzscore profile features "+str(perP)
    Wtrain, Sigt, Ht = np.linalg.svd(zscores.T, full_matrices=False)
    
    dimP = min(len(Sigt), len(np.where(np.cumsum(Sigt)/np.sum(Sigt) < perP)[0])+1)
    print "new dim", dimP

    Wtrain = Wtrain[:, :dimP]
    Sigt = Sigt[:dimP]
    Htrain = Ht[: dimP, :]
    
    zscores = Wtrain.T
    
    
elif "--NMF" in sys.argv:
    # use NMF when less features than data?
    print "NMF..."
    perP = sys.argv[sys.argv.index('--NMF')+1]
    if perP.isdigit():
        dimP = int(perP)
    else:
        dimP = len(Ptrain)
    from sklearn.decomposition import NMF
    model = NMF(n_components=dimP)
    model.fit(zscores)
    Htrain = model.components_
    Wtrain = model.transform(zscores)
    
if "--outname" in sys.argv:
    outname = sys.argv[sys.argv.index("--outname")+1]
else:
    print "Please provide outname"
    sys.exit()
    
corrmat = np.ones((len(motnames),len(motnames)))


if '--determinetop' in sys.argv:
    sigcut = float(sys.argv[sys.argv.index('--determinetop')+1])
    sigcut = int(float(np.shape(zscores)[0]*np.shape(zscores)[1])*sigcut)
    cutvalue = np.sort(zscores.flatten())[-sigcut]
    print cutvalue
    topvalues = []
    for zscore in zscores.T:
        if len(np.nonzero(zscore > cutvalue)[0]) > 0:
            topvalues.append(zscore > cutvalue)
        else:
            topvalues.append(zscore == np.amax(zscore))
            #print np.nonzero(topvalues[-1])[0]
        
elif '--topnum' in sys.argv:
    topnum = int(sys.argv[sys.argv.index('--topnum')+1])
    topvalues = []
    for zscore in zscores.T:
        topvalues.append(zscore > -np.sort(-zscore)[topnum])

elif '--topfraction' in sys.argv:
    topfraction = float(sys.argv[sys.argv.index('--topfraction')+1])
    topnum = int(float(np.shape(zscores)[0]) * topfraction)
    topvalues = []
    for zscore in zscores.T:
        topvalues.append(zscore > -np.sort(-zscore)[topnum])

                      
elif '--topcut' in sys.argv:
    cutvalue = float(sys.argv[sys.argv.index('--topcut')+1])
    topvalues = []
    for zscore in zscores.T:
        if len(np.nonzero(zscore > cutvalue)[0]) > 0:
            topvalues.append(zscore > cutvalue)
        else:
            topvalues.append(zscore == np.amax(zscore))

else:
    topvalues = np.ones(np.shape(zscores.T), dtype = int) > 0




if "--pearson" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cout = pearsonr(zvec1[maskv], zvec2[maskv])[0]
        return cout

elif '--eigenvalue' in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        ab = np.sum(zvec1[maskv]* zvec2[maskv])
        aa = np.sum(zvec1[maskv]* zvec1[maskv])
        bb = np.sum(zvec2[maskv]* zvec2[maskv])
        cout = ((aa+bb)/2.+np.sqrt((aa-bb)**2/4.+ab**2))/((aa+bb)/2.-np.sqrt((aa-bb)**2/4.+ab**2))
        return cout

elif '--eigenvector' in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cov = np.zeros((2,2))
        cov[0,1] = cov[1,0] = np.sum(zvec1[maskv]* zvec2[maskv])
        cov[0,0] = np.sum(zvec1[maskv]* zvec1[maskv])
        cov[1,1] = np.sum(zvec2[maskv]* zvec2[maskv])
        csig, cevec = np.linalg.eig(cov)
        return np.sum(cevec[:,0]*np.array([1,1])/np.sqrt(2.))

elif "--pearsonweighted" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        weight = np.zeros(len(zvec1))
        weight[maskv1] += 0.5
        weight[maskv2] += 0.5
        maskv = maskv1 | maskv2
        cout = np.sum((zvec1[maskv]-np.mean(zvec1[maskv]))*weight[maskv]*(zvec2[maskv]-np.mean(zvec2[maskv])))/(np.sum(weight[maskv])*np.std(zvec1[maskv])*np.std(zvec2[maskv]))
        return cout

elif '--normdot' in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cout = np.dot(zvec1[maskv], zvec2[maskv])**2/(np.dot(zvec1[maskv], zvec1[maskv])*np.dot(zvec2[maskv], zvec2[maskv]))
        return cout

elif "--dotproduct" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cout = np.dot(zvec1[maskv], zvec2[maskv])
        return cout

elif '--mdot' in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cout = 2.*np.dot(zvec1[maskv], zvec2[maskv])/(np.dot(zvec1[maskv], zvec1[maskv])+np.dot(zvec2[maskv], zvec2[maskv]))
        return cout
    
elif "--spearman" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cout = spearmanr(zvec1[maskv], zvec2[maskv])[0]
        return cout
    
elif "--cosine" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cout = cosine(zvec1[maskv], zvec2[maskv])
        return cout
    
elif "--coefofdet" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        maskv = maskv1 | maskv2
        cout = r2_score(zvec1[maskv], zvec2[maskv])
        return cout

elif '--meansquaredeviation' in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        # top deviation devided by self deviation
        top1 = maskv1
        top2 = maskv2
        cout = 1.- (np.mean((zvec1[top1]-zvec2[top1])**2)+np.mean((zvec1[top2]-zvec2[top2])**2))/(np.mean((zvec1[top1]-np.mean(zvec1[top1])**2))+np.mean((zvec2[top2]-np.mean(zvec2[top2])**2)))
        return cout

elif '--meanabsdeviation' in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        top1 = maskv1
        top2 = maskv2
        cout = 1.-(np.mean(np.absolute(zvec1[top1]-zvec2[top1]))+np.mean(np.absolute(zvec1[top2]-zvec2[top2])))/(np.mean(np.absolute(zvec1[top1]-np.mean(zvec1[top1])))+np.mean(np.absolute(zvec2[top2]-np.mean(zvec2[top2]))))
        return cout


###### topset and classes go wrong somehow!!!!! , maybe because of topcut5
elif "--topset" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        #print len(np.intersect1d(np.arange(len(zvec1))[maskv1], np.arange(len(zvec1))[maskv2])), np.sum(maskv1)
        cout = float(len(np.intersect1d(np.arange(len(zvec1))[maskv1], np.arange(len(zvec1))[maskv2])))/np.sum(maskv1)
        #print cout
        return cout    




elif "--topsetweighted" in sys.argv:
    def cmeth(zvec1, zvec2, maskv1, maskv2):
        setAen = np.arange(len(zvec1))[maskv1]
        setBen = np.arange(len(zvec2))[maskv2]
        overlapC = np.intersect1d(setAen, setBen)
        unionE = np.union1d(setAen, setBen)
        cout = float(np.sum(zvec1[overlapC]+zvec2[overlap])/2.)/float(np.sum(zvec1[unionE]+zvec2[unionE])/2.)
        return cout  
  
else: 
    print "Please provide correlation method"
    sys.exit()


    
for m in range(len(motnames)):
    #print m
    for n in range(m+1,len(motnames)):
        corrmat[m,n] = corrmat[n,m] = cmeth(zscores[:,m], zscores[:,n], topvalues[m], topvalues[n])
        #print corrmat[m,n]
        #if corrmat[m,n] > 2.:
            #print corrmat[m,n] 
corrmat = np.nan_to_num(corrmat)
if '--savetxt' in sys.argv:
    np.savetxt(outname, corrmat, fmt = '%2.4f', header = ' '.join(motnames))

if "--print_top" in sys.argv:
    top = int(sys.argv[sys.argv.index("--print_top")+1])
    corrmatsort = np.argsort(corrmat)
    for i, sort in enumerate(corrmatsort):
        print motnames[i],'|', ' '.join(motnames[sort[-top-1:-1]]), corrmat[i,sort[-top-1:-1]]
        
    
    
