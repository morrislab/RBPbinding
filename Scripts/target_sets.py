import numpy as np
import sys, os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.cluster import hierarchy
import logomaker as lm
from scipy.stats import skew
from matplotlib.colors import ListedColormap
from matplotlib import cm
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import pearsonr

def ic(pw):
    icout = pw*np.log2(pw/0.25)
    icout[icout < -2] = -2
    return icout

def plotpwm(pwm, ax):
    #print pwm
    pwm = pd.DataFrame({'A':pwm[:,0],'C':pwm[:,1], 'G':pwm[:,2], 'U':pwm[:,3]})        
    pwm = ic(pwm)
    lm.Logo(pwm, ax = ax)
    #ax.set_yticks([0,1,2])
    #ax.set_yticklabels([0,1,2])
    ax.set_ylim([0,2])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params('both', left=False, bottom = False, labelleft = False, labelbottom = False)
    return ax


def readinpwm(pwfile):
    obj = open(pwfile,'r').readlines()
    pws = []
    pwname = []
    for l, line in enumerate(obj):
        line = line.strip().split()
        if len(line) != 0:
            if line[0] == 'Motif':
                if len(pwname) > 0:
                    pws.append(np.array(pw, dtype = float))
                pwname.append(line[1])
                pw = []
            if line[0].isdigit():
                pw.append(line[1:])
    pws.append(np.array(pw, dtype = float))
    return pwname, pws





efile = np.load(sys.argv[1])
correlationmat = efile['correlationmat']
bindingmat = efile['bindingmat']
rbpnames=efile['rbpnames']
genes=efile['genes']

pset = np.genfromtxt(sys.argv[2], dtype = str)
outname = os.path.splitext(sys.argv[2])[0]
crow, srow = sys.argv[3].split(',')
pnames = np.genfromtxt(sys.argv[4], dtype = str)
pwmfile = sys.argv[5]
outname += '_pos'+srow

folder = sys.argv[2].rsplit('.',1)[0]

pwmname, pwms = readinpwm(pwmfile)

pset = pset[pset[:, int(srow)] == 'True']

colors = np.array(['grey', 'coral', 'cornflowerblue'])
rbpcolors = np.zeros(len(pset), dtype = int)
rbpcolors[pset[:, int(crow)].astype(float) > 0] = 1
rbpcolors[pset[:, int(crow)].astype(float) < 0] = 2
rpbcolors = colors[rbpcolors]
zscores = pset[:, int(crow)]
pset = pset[:, 0]

zsort = np.argsort(zscores.astype(float))

pset = pset[zsort]
zscores = zscores[zsort]
rbpcolors = rbpcolors[zsort]


nameset = []
setprotnames = []
setgenenames = []
for p, pr in enumerate(pset):
    if '--double' in sys.argv:
        setprotnames.append(pnames[list(pnames[:, 0]).index(pr.split('--')[0]), 3]+'-'+pnames[list(pnames[:, 0]).index(pr.split('--')[1]), 3])
        setgenenames.append(pnames[list(pnames[:, 0]).index(pr.split('--')[0]), 2]+'-'+pnames[list(pnames[:, 0]).index(pr.split('--')[1]), 2])
    else:
        setprotnames.append(pnames[list(pnames[:, 0]).index(pr), 3])
        setgenenames.append(pnames[list(pnames[:, 0]).index(pr), 2].split('.')[0])
    nameset.append(setprotnames[-1])

bexpsort = []
pwsort = []
for gname in pset:
    bexpsort.append(list(rbpnames).index(gname))
    if '--double' in sys.argv:
        prs = gname.split('--')
        pwsort.append(np.concatenate([pwms[list(pwmname).index(prs[0])], np.ones((2,4))*0.25, pwms[list(pwmname).index(prs[1])]],axis =0))
    else:
        pwsort.append(pwms[list(pwmname).index(gname)])
    
pwms = pwsort
bindingmat = bindingmat[bexpsort]
correlationmat = correlationmat[bexpsort]

targetgenes = []
nontargetgenes = []
genecorrelation = []
targetscore = []

if '--cut' in sys.argv:
    ccut = float(sys.argv[sys.argv.index('--cut')+1])
    folder += '_cut'+str(ccut)
    for p, prot in enumerate(pset):
        targ = (np.sign(float(zscores[p]))*correlationmat[p]>=ccut)*bindingmat[p]
        genecorrelation.append(np.where(targ)[0])
        targetgenes.append(genes[targ])
        targetscore.append(np.sign(float(zscores[p]))*correlationmat[p][targ])
        nontargetgenes.append(genes[~targ])
        
elif '--max_dist' in sys.argv:
    folder += '_maxdist'
    for p, prot in enumerate(pset):
        cormat = np.sign(float(zscores[p]))*correlationmat[p]
        sortc = np.argsort(cormat)
        bmat = bindingmat[p][sortc]
        
        cdfmat = np.zeros((2, len(bmat)))
        cdfmat[0,bmat] = np.ones(int(np.sum(bmat)))/np.sum(bmat)
        cdfmat[1,~bmat] = np.ones(int(np.sum(~bmat)))/np.sum(~bmat)
        cdfmat = np.cumsum(cdfmat, axis = 1)
        
        cdfdist = np.diff(cdfmat, axis = 0)[0]
        plt.plot(cormat[sortc], cdfmat[0])
        plt.plot(cormat[sortc], cdfmat[1])
        plt.plot(cormat[sortc], cdfdist)
        plt.show()
        cdfdist[cormat[sortc] <= 0] = 0
        cormax = cormat[sortc][np.argmax(cdfdist)]
        print setprotnames[p], cormax
        genecorrelation.append(np.where((cormat>=cormax) * bindingmat[p])[0])
        targetscore.append(cormat[genecorrelation[-1]])
        targetgenes.append(genes[genecorrelation[-1]])
        nontargetgenes.append(genes[np.delete(np.arange(len(genes)),genecorrelation[-1])])
    
   
np.savez_compressed(folder+'-targetgenes.npz', pset = pset, genenames = setgenenames, protnames = setprotnames, targets = targetgenes, targetscores = targetscore, nontargetgenes= nontargetgenes, zscores = np.array(zscores, dtype = float))

if '--plot' in sys.argv:
    fractionmat = np.zeros((len(nameset), len(nameset)))
    overlapmat = np.zeros((len(nameset), len(nameset)))
    clustermat = np.zeros((len(nameset), len(nameset)))
    for g, gc in enumerate(genecorrelation):
        for h in range(len(genecorrelation)):
            overlapmat[g, h] = len(np.intersect1d(gc, genecorrelation[h]))
            fractionmat[g, h] = len(np.intersect1d(gc, genecorrelation[h]))/float(len(gc))
            clustermat[g,h] = len(np.intersect1d(gc, genecorrelation[h]))/float(len(np.union1d(gc, genecorrelation[h])))

 
    cmapgs = np.concatenate([[[1.,1.,1.,1.]] , [cm.Purples(0.1) for i in range(9)], [cm.Purples(0.4) for i in range(10)], [cm.Purples(0.55) for i in range(6)] , [cm.Purples(0.7) for i in range(4)], [cm.Purples(1.) for i in range(10)]], axis = 0)

    cmapgs2 = np.concatenate([[[1.,1.,1.,1.] for i in range(1)] , [cm.Purples(0.1) for i in range(9)], [cm.Purples(0.3) for i in range(10)],[cm.Purples(0.5) for i in range(10)],  [cm.Purples(0.7) for i in range(10)], [cm.Purples(1.) for i in range(10)]], axis = 0)

    colmap = ListedColormap(cmapgs)
    colmap2 = ListedColormap(cmapgs2)
    tcolmap = ListedColormap(colors)



    fig = plt.figure(figsize=(10,8), dpi = 100)

    aden = fig.add_subplot(261)
    aden.set_position([0.1, 0.1, 0.1 , 0.8])
    aden.spines['top'].set_visible(False)
    aden.spines['bottom'].set_visible(False)
    aden.spines['left'].set_visible(False)
    aden.spines['right'].set_visible(False)
    aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Z = linkage(1.-clustermat, 'single')
    dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
    sorting = dn['leaves']


    ac = fig.add_subplot(266)
    ac.set_position([0.7, 0.925, 0.15 ,0.05])
    ac.imshow([np.concatenate([np.linspace(0,2, 20), np.linspace(2,2.6, 10), np.linspace(2.6,3, 10), np.linspace(3,4, 10)])], origin = 'lower', cmap = colmap, aspect = 'auto')
    ac.set_xticks([0,9,19,29,39,49])
    ac.set_xticklabels([0,10, 100, 500, 1000, 10000], rotation = 60)
    ac.spines['top'].set_visible(False)
    ac.spines['left'].set_visible(False)
    ac.spines['right'].set_visible(False)
    ac.spines['bottom'].set_visible(False)        
    ac.tick_params(left = False, labelleft = False, right = False, labelright = False, bottom = False, labelbottom = False, top = True, labeltop = True)
    #ac.set_title('Target overlap')

    pheight = 0.8/float(len(nameset))
    for n, s in enumerate(sorting):
        print n
        axp = fig.add_subplot(len(nameset), 1, n+1)
        axp.set_position([0.2, 0.1+pheight*(n+0.1), 0.075, pheight*0.8])
        plotpwm(pwms[s], axp)
        

    # imshow for rbpcolors
    axt = fig.add_subplot(265)
    axt.set_position([0.275, 0.1, 0.05, 0.8])
    axt.spines['top'].set_visible(False)
    axt.spines['bottom'].set_visible(False)
    axt.spines['left'].set_visible(False)
    axt.spines['right'].set_visible(False)
    axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False, labelbottom = False, labeltop = True)
    imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
    axt.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
    axt.grid(color='silver', linewidth=0.5, which = 'minor')
    axt.set_xticks([-0.1])
    axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'bottom', ha = 'left')
    for i in range(len(sorting)):
        axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center')


    ax = fig.add_subplot(263)
    ax.set_position([0.325, 0.1, 0.575, 0.8])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(which = 'both', labelleft = False, left = False, labelbottom = False, bottom = False, labelright = True)
    im = ax.imshow(np.log10(overlapmat+1)[sorting][:, sorting], origin = 'lower', cmap = colmap, aspect = 'auto', vmin = 0, vmax = 4)
    ax.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
    ax.set_xticks( np.arange(len(nameset)+1)-0.5, minor = True)
    ax.grid(color='silver', linewidth=0.5, which = 'minor')
    ax.set_yticks(np.arange(len(nameset)))
    ax.set_yticklabels(np.array(nameset)[sorting]) #rotation = 60, va = 'top', ha = 'right') 

    if '--plot' in sys.argv:
        fig.savefig(folder+'-targetoverlap.jpg', bbox_inches = 'tight', dpi = 300)
        plt.show()














    fig = plt.figure(figsize=(10,8), dpi = 100)

    aden = fig.add_subplot(261)
    aden.set_position([0.1, 0.1, 0.1 , 0.8])
    aden.spines['top'].set_visible(False)
    aden.spines['bottom'].set_visible(False)
    aden.spines['left'].set_visible(False)
    aden.spines['right'].set_visible(False)
    aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Z = linkage(1.-clustermat, 'single')
    dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
    sorting = dn['leaves']


    ac = fig.add_subplot(266)
    ac.set_position([0.7, 0.925, 0.15 ,0.05])
    ac.imshow([np.linspace(0,1, 50)], origin = 'lower', cmap = colmap2, aspect = 'auto')
    ac.set_xticks([0,24,49])
    ac.set_xticklabels([0,0.5,1.], rotation = 60)
    ac.spines['top'].set_visible(False)
    ac.spines['left'].set_visible(False)
    ac.spines['right'].set_visible(False)
    ac.spines['bottom'].set_visible(False)        
    ac.tick_params(left = False, labelleft = False, right = False, labelright = False, bottom = False, labelbottom = False, top = True, labeltop = True)
    #ac.set_title('Target overlap')

    pheight = 0.8/float(len(nameset))
    for n, s in enumerate(sorting):
        print n
        axp = fig.add_subplot(len(nameset), 1, n+1)
        axp.set_position([0.2, 0.1+pheight*(n+0.1), 0.075, pheight*0.8])
        plotpwm(pwms[s], axp)
        

    # imshow for rbpcolors
    axt = fig.add_subplot(265)
    axt.set_position([0.275, 0.1, 0.05, 0.8])
    axt.spines['top'].set_visible(False)
    axt.spines['bottom'].set_visible(False)
    axt.spines['left'].set_visible(False)
    axt.spines['right'].set_visible(False)
    axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False, labelbottom = False, labeltop = True)
    imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
    axt.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
    axt.grid(color='silver', linewidth=0.5, which = 'minor')
    axt.set_xticks([-0.1])
    axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'bottom', ha = 'left')
    for i in range(len(sorting)):
        axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center')


    ax = fig.add_subplot(263)
    ax.set_position([0.325, 0.1, 0.475, 0.8])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(which = 'both', labelleft = False, left = False, labelbottom = False, bottom = False, labelright = True)
    im = ax.imshow(fractionmat[sorting][:, sorting], origin = 'lower', cmap = colmap2, aspect = 'auto', vmin = 0, vmax = 1)
    ax.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
    ax.set_xticks( np.arange(len(nameset)+1)-0.5, minor = True)
    ax.grid(color='silver', linewidth=0.5, which = 'minor')
    ylim = ax.get_ylim()

    axbar = fig.add_subplot(267)
    axbar.set_position([0.8,0.1,0.1,0.8])

    print np.diagonal(overlapmat)
    axbar.barh(np.arange(len(nameset)), np.diagonal(overlapmat)[sorting], height = 0.75, color = 'darkslateblue', edgecolor = 'k') 
    axbar.spines['right'].set_visible(False)
    axbar.spines['top'].set_visible(False)
    axbar.spines['left'].set_visible(False)
    axbar.set_xlabel('Number\ntargets')
    axbar.tick_params(left = False, labelleft = False, labelright = True)
    axbar.set_yticks(np.arange(len(nameset)))
    axbar.set_yticklabels(np.array(nameset)[sorting])
    axbar.set_ylim(ylim)


    if '--plot' in sys.argv:
        fig.savefig(folder+'-targetfractionoverlap.jpg', bbox_inches ='tight', dpi =300)
        plt.show()



