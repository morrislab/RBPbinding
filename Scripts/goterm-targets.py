import numpy as np
import sys, os
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.cluster import hierarchy
import logomaker as lm
from scipy.stats import skew, fisher_exact
from matplotlib.colors import ListedColormap
from matplotlib import cm
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

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


pset = np.load(sys.argv[1], allow_pickle = True)
#print pset.files
setgenenames = pset['genenames']
setprotnames = pset['protnames']
targetsets = pset['targets']
targetscores = pset['targetscores']
nontargetsets = pset['nontargetgenes']
zscores = pset['zscores']
pset = pset['pset']

folder = os.path.splitext(sys.argv[1])[0]

pwmfile = sys.argv[2]

pwmname, pwms = readinpwm(pwmfile)

colors = np.array(['grey', 'coral', 'cornflowerblue'])
rbpcolors = np.zeros(len(pset), dtype = int)
rbpcolors[zscores > 0] = 1
rbpcolors[zscores < 0] = 2


nameset = []
for p, pr in enumerate(setprotnames):
    if '--double' in sys.argv:
        nameset.append(pr)
    else:
        nameset.append(pr+' ['+setgenenames[p]+']')
    
pwsort = []
for gname in pset:
    if '--double' in sys.argv:
        pwsort.append(np.concatenate([pwms[list(pwmname).index(gname.split('--')[0])], np.ones((2,4))*0.25, pwms[list(pwmname).index(gname.split('--')[1])]],axis =0))
    else:
        pwsort.append(pwms[list(pwmname).index(gname)])
pwms = pwsort

gofile = np.load(sys.argv[3])
genes = gofile['genes']

gopmat = gofile['gopmat']
gofmat = gofile['gofmat']
gocmat = gofile['gocmat']

gotermp = gofile['gotermp']
goidp = gofile['goidp']

gotermf = gofile['gotermf']
goidf = gofile['goidf']

gotermc = gofile['gotermc']
goidc = gofile['goidc']


sortp = np.argsort(gotermp)
sortc = np.argsort(gotermc)
sortf = np.argsort(gotermf)

gotermp = gotermp[sortp]
gotermf = gotermf[sortf]
gotermc = gotermc[sortc]

gopmat = gopmat[:, sortp]
gocmat = gocmat[:, sortc]
gofmat = gofmat[:, sortf]



def enrichment(tid, ntid, mat):
    if len(ntid) == 0:
        ntid = np.delete(np.arange(len(mat)), tid)
    tsetlen = float(len(tid))
    nsetlen = float(len(ntid))
    
    tgset = np.sum(mat[tid], axis = 0)
    tngset = np.sum(mat[ntid], axis = 0)
    
    tgsetmean = tgset/tsetlen
    tngsetmean = tngset/nsetlen
    
    #tgsetstd = np.sqrt(tgsetmean-2.*tgsetmean**2+tgsetmean**3)
    tngsetstd = np.sqrt(tngsetmean-2.*tngsetmean**2+tngsetmean**3)
    
    hasobservation = np.where((tgset > 1))[0]
    consider = np.where((tgset[hasobservation] > 1)*(tgsetmean[hasobservation] > tngsetmean[hasobservation]+tngsetstd[hasobservation]))[0]
    
    if len(consider) > 0:
        pvals = np.ones(len(hasobservation))
        for c in consider:
            odr, pval = fisher_exact([[tgset[hasobservation[c]],tngset[hasobservation[c]]],[tsetlen-tgset[hasobservation[c]],nsetlen-tngset[hasobservation[c]]]], 'greater')
            pvals[c] = pval
        if '--singletesting' in sys.argv:
            msig, multip = pvals < 0.05, pvals
        else:
            msig, multip, carr, carr2 = multipletests(pvals, alpha = 0.05, method = 'fdr_bh')
        #print multip[multip< 1], pvals[multip< 1]
        return hasobservation[msig], tgsetmean[hasobservation[msig]]
    else:
        return [], []


sigtermp = []
sigtermf = []
sigtermc = []

sigtermpn = []
sigtermfn = []
sigtermcn = []
for p, ps in enumerate(pset):
    tid = np.where(np.isin(genes,targetsets[p]))[0]
    ntid = nontargetsets[p]
    if len(ntid) > 0:
        ntid = np.where(np.isin(genes,ntid))[0]
    if '--print_individual' in sys.argv:
        print '\n',setprotnames[p], len(tid)
        print genes[tid]
        numocc = np.sum(gopmat[tid], axis = 0)
        print 'P'
        for s, a in enumerate(gotermp[numocc >=1]):
            print a, numocc[numocc>=1][s]
        numocc = np.sum(gocmat[tid], axis = 0)
        print 'C'
        for s, a in enumerate(gotermc[numocc >=1]):
            print a, numocc[numocc>=1][s]
        numocc = np.sum(gofmat[tid], axis = 0)
        print 'F'
        for s, a in enumerate(gotermf[numocc >=1]):
            print a, numocc[numocc>=1][s]
        print '\n\n\n'
    enrichedp, fracp = enrichment(tid, ntid, gopmat)
    enrichedc, fracc = enrichment(tid, ntid, gocmat)
    enrichedf, fracf = enrichment(tid, ntid, gofmat)
   
    sigtermp.append(gotermp[enrichedp])
    sigtermc.append(gotermc[enrichedc])
    sigtermf.append(gotermf[enrichedf])
    
    if (len(sigtermp[-1]) > 0) or (len(sigtermf[-1]) > 0) or (len(sigtermc[-1]) > 0):
        print '\n', nameset[p]
    
    if len(sigtermp[-1]) > 0:
        print 'P'
        print '\n'.join(sigtermp[-1])
        
    if len(sigtermc[-1]) > 0:
        print 'C'
        print '\n'.join(sigtermc[-1])
        
    if len(sigtermf[-1]) > 0:
        print 'F'
        print '\n'.join(sigtermf[-1])
        
    
    sigtermpn.append(fracp)
    sigtermcn.append(fracc)
    sigtermfn.append(fracf)


unsigp = np.unique(np.concatenate(sigtermp))
unsigf = np.unique(np.concatenate(sigtermf))
unsigc = np.unique(np.concatenate(sigtermc))

#print unsigp, unsigf, unsigc

sigpmat = np.zeros((len(pset), len(unsigp))) - 0.02
sigfmat = np.zeros((len(pset), len(unsigf))) - 0.02
sigcmat = np.zeros((len(pset), len(unsigc))) - 0.02

for p in range(len(pset)):
    sigpmat[p, np.isin(unsigp, sigtermp[p])] = sigtermpn[p]
    sigfmat[p, np.isin(unsigf, sigtermf[p])] = sigtermfn[p]
    sigcmat[p, np.isin(unsigc, sigtermc[p])] = sigtermcn[p]


print np.shape(sigpmat), np.shape(sigfmat), np.shape(sigcmat)

tcolmap = ListedColormap(colors)

if np.shape(sigpmat)[1] > 1:
    pmask = np.where(np.sum(sigpmat > -0.02, axis = 1) > 0)[0]

    fig = plt.figure(figsize=(len(sigpmat[0])*0.3+2.5,1+len(pmask)*0.3), dpi = 100)

    dw = 0.8/(len(sigpmat[0])*0.3 +2.5)
    aden = fig.add_subplot(261)
    aden.set_position([0.1, 0.1, dw*1. , 0.8])
    aden.spines['top'].set_visible(False)
    aden.spines['bottom'].set_visible(False)
    aden.spines['left'].set_visible(False)
    aden.spines['right'].set_visible(False)
    aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Z = linkage(cdist(sigpmat[pmask], sigpmat[pmask], 'euclidean'), 'single')
    dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
    sorting = dn['leaves']
    sorting = pmask[sorting]

    adent = fig.add_subplot(262)
    adent.set_position([0.1+2.5*dw, 0.9, 0.8-2.5*dw , 0.075])
    adent.spines['top'].set_visible(False)
    adent.spines['bottom'].set_visible(False)
    adent.spines['left'].set_visible(False)
    adent.spines['right'].set_visible(False)
    adent.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Zg = linkage(cdist(sigpmat[pmask].T, sigpmat[pmask].T, 'euclidean'), 'single')
    dng = dendrogram(Zg, orientation = 'top', color_threshold=0, above_threshold_color='k', ax = adent)
    gsorting = dng['leaves']


    pheight = 0.8/float(len(sorting))
    for n, s in enumerate(sorting):
        print n
        axp = fig.add_subplot(len(sorting), 1, n+1)
        axp.set_position([0.1+1.*dw, 0.1+pheight*(n+0.1), 1.*dw, pheight*0.8])
        plotpwm(pwms[s], axp)
        

    # imshow for rbpcolors
    axt = fig.add_subplot(265)
    axt.set_position([0.1+2*dw, 0.1, 0.5*dw, 0.8])
    axt.spines['top'].set_visible(False)
    axt.spines['bottom'].set_visible(False)
    axt.spines['left'].set_visible(False)
    axt.spines['right'].set_visible(False)
    axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False, labelbottom = True, labeltop = False)
    imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
    axt.set_yticks(np.arange(len(sorting)+1)-0.5, minor = True)
    axt.grid(color='silver', linewidth=0.5, which = 'minor')
    axt.set_xticks([0.1])
    axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'top', ha = 'right')
    for i in range(len(sorting)):
        axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center')


    ax = fig.add_subplot(263)
    ax.set_position([0.1+2.5*dw, 0.1, 0.8-2.5*dw, 0.8])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(which = 'both', labelleft = False, left = False, labelright = True, bottom = False)
    im = ax.imshow(sigpmat[sorting][:, gsorting], origin = 'lower', cmap = 'Greys', aspect = 'auto')
    ax.set_yticks(np.arange(len(sorting)+1)-0.5, minor = True)
    ax.set_xticks( np.arange(len(unsigp)+1)-0.5, minor = True)
    ax.grid(color='silver', linewidth=0.5, which = 'minor')
    ax.set_yticks(np.arange(len(sorting)))
    ax.set_yticklabels(np.array(nameset)[sorting])

    ax.set_xticks(np.arange(len(unsigp))+0.25)
    ax.set_xticklabels(np.array(unsigp)[gsorting], va = 'top', ha = 'right', rotation = 60)

    sortedmat = sigpmat[sorting][:,gsorting]
    for i in range(len(sortedmat)):
        for j in range(len(sortedmat[0])):
            if sortedmat[i,j] > 0:
                ax.text(j,i, str(np.around(sortedmat[i,j]*100,1)), color = 'white', va = 'center', ha = 'center')


    if '--savefig' in sys.argv:
        fig.savefig(os.path.splitext(sys.argv[1])[0]+'-goprocess-targets.jpg', bbox_inches = 'tight', dpi = 300)
        print os.path.splitext(sys.argv[1])[0]+'-goprocess-targets.jpg' 

if np.shape(sigfmat)[1] > 1:
    fmask = np.where(np.sum(sigfmat > -0.02, axis = 1) > 0)[0]
    print fmask
    fig2 = plt.figure(figsize=(len(sigfmat[0])*0.3+2.5,len(fmask)*0.3+1), dpi = 100)
    dw = 0.8/(len(sigfmat[0])*0.3 +2.5)

    aden = fig2.add_subplot(261)
    aden.set_position([0.1, 0.1, 1.*dw , 0.8])
    aden.spines['top'].set_visible(False)
    aden.spines['bottom'].set_visible(False)
    aden.spines['left'].set_visible(False)
    aden.spines['right'].set_visible(False)
    aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Z = linkage(cdist(sigfmat[fmask], sigfmat[fmask], 'euclidean'), 'single')
    dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
    sorting = dn['leaves']
    sorting = fmask[sorting]

    adent = fig2.add_subplot(262)
    adent.set_position([0.1 + 2.5*dw, 0.9, 0.8-2.5*dw , 0.075])
    adent.spines['top'].set_visible(False)
    adent.spines['bottom'].set_visible(False)
    adent.spines['left'].set_visible(False)
    adent.spines['right'].set_visible(False)
    adent.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Zg = linkage(cdist(sigfmat[fmask].T, sigfmat[fmask].T, 'euclidean'), 'single')
    dng = dendrogram(Zg, orientation = 'top', color_threshold=0, above_threshold_color='k', ax = adent)
    gsorting = dng['leaves']


    pheight = 0.8/float(len(sorting))
    for n, s in enumerate(sorting):
        print n
        axp = fig2.add_subplot(len(sorting), 1, n+1)
        axp.set_position([0.1+1.*dw, 0.1+pheight*(n+0.1), 1.*dw, pheight*0.8])
        plotpwm(pwms[s], axp)
        

    # imshow for rbpcolors
    axt = fig2.add_subplot(265)
    axt.set_position([0.1 + 2.*dw, 0.1, 0.5*dw, 0.8])
    axt.spines['top'].set_visible(False)
    axt.spines['bottom'].set_visible(False)
    axt.spines['left'].set_visible(False)
    axt.spines['right'].set_visible(False)
    axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False, labelbottom = True, labeltop = False)
    imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
    axt.set_yticks(np.arange(len(sorting)+1)-0.5, minor = True)
    axt.grid(color='silver', linewidth=0.5, which = 'minor')
    axt.set_xticks([0.1])
    axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'top', ha = 'right')
    for i in range(len(sorting)):
        axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center')


    ax = fig2.add_subplot(263)
    ax.set_position([0.1 + 2.5*dw, 0.1, 0.8-2.5*dw, 0.8])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(which = 'both', labelleft = False, left = False, labelright = True, bottom = False)
    im = ax.imshow(sigfmat[sorting][:, gsorting], origin = 'lower', cmap = 'Greys', aspect = 'auto')
    ax.set_yticks(np.arange(len(sorting)+1)-0.5, minor = True)
    ax.set_xticks( np.arange(len(unsigf)+1)-0.5, minor = True)
    ax.grid(color='silver', linewidth=0.5, which = 'minor')
    ax.set_yticks(np.arange(len(sorting)))
    ax.set_yticklabels(np.array(nameset)[sorting])

    ax.set_xticks(np.arange(len(unsigf))+0.25)
    ax.set_xticklabels(np.array(unsigf)[gsorting], va = 'top', ha = 'right', rotation = 60)

    sortedmat = sigfmat[sorting][:,gsorting]
    for i in range(len(sortedmat)):
        for j in range(len(sortedmat[0])):
            if sortedmat[i,j] > 0:
                ax.text(j,i, str(np.around(sortedmat[i,j]*100,1)), color = 'white', va = 'center', ha = 'center')


    if '--savefig' in sys.argv:
        fig2.savefig(os.path.splitext(sys.argv[1])[0]+'-gofunction-targets.jpg', bbox_inches = 'tight', dpi = 300)







if np.shape(sigcmat)[1] > 1:
    cmask = np.where(np.sum(sigcmat > -0.02, axis = 1) > 0)[0]
    fig3 = plt.figure(figsize=(len(sigcmat[0])*0.5,len(cmask)*0.5), dpi = 100)


    aden = fig3.add_subplot(261)
    aden.set_position([0.1, 0.1, 0.1 , 0.8])
    aden.spines['top'].set_visible(False)
    aden.spines['bottom'].set_visible(False)
    aden.spines['left'].set_visible(False)
    aden.spines['right'].set_visible(False)
    aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Z = linkage(cdist(sigcmat[cmask], sigcmat[cmask], 'euclidean'), 'single')
    dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
    sorting = dn['leaves']
    sorting = cmask[sorting]

    adent = fig3.add_subplot(262)
    adent.set_position([0.325, 0.9, 0.575 , 0.075])
    adent.spines['top'].set_visible(False)
    adent.spines['bottom'].set_visible(False)
    adent.spines['left'].set_visible(False)
    adent.spines['right'].set_visible(False)
    adent.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
    Zg = linkage(cdist(sigcmat[cmask].T, sigcmat[cmask].T, 'euclidean'), 'single')
    dng = dendrogram(Zg, orientation = 'top', color_threshold=0, above_threshold_color='k', ax = adent)
    gsorting = dng['leaves']


    pheight = 0.8/float(len(sorting))
    for n, s in enumerate(sorting):
        print n
        axp = fig3.add_subplot(len(sorting), 1, n+1)
        axp.set_position([0.2, 0.1+pheight*(n+0.1), 0.075, pheight*0.8])
        plotpwm(pwms[s], axp)
        

    # imshow for rbpcolors
    axt = fig3.add_subplot(265)
    axt.set_position([0.275, 0.1, 0.05, 0.8])
    axt.spines['top'].set_visible(False)
    axt.spines['bottom'].set_visible(False)
    axt.spines['left'].set_visible(False)
    axt.spines['right'].set_visible(False)
    axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False, labelbottom = True, labeltop = False)
    imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
    axt.set_yticks(np.arange(len(sorting)+1)-0.5, minor = True)
    axt.grid(color='silver', linewidth=0.5, which = 'minor')
    axt.set_xticks([0.1])
    axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'top', ha = 'right')
    for i in range(len(sorting)):
        axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center')


    ax = fig3.add_subplot(263)
    ax.set_position([0.325, 0.1, 0.575, 0.8])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(which = 'both', labelleft = False, left = False, labelright = True, bottom = False)
    im = ax.imshow(sigcmat[sorting][:, gsorting], origin = 'lower', cmap = 'Greys', aspect = 'auto')
    ax.set_yticks(np.arange(len(sorting)+1)-0.5, minor = True)
    ax.set_xticks( np.arange(len(unsigc)+1)-0.5, minor = True)
    ax.grid(color='silver', linewidth=0.5, which = 'minor')
    ax.set_yticks(np.arange(len(sorting)))
    ax.set_yticklabels(np.array(nameset)[sorting])

    ax.set_xticks(np.arange(len(unsigc))+0.25)
    ax.set_xticklabels(np.array(unsigc)[gsorting], va = 'top', ha = 'right', rotation = 60)

    sortedmat = sigcmat[sorting][:,gsorting]
    for i in range(len(sortedmat)):
        for j in range(len(sortedmat[0])):
            if sortedmat[i,j] > 0:
                ax.text(j,i, str(np.around(sortedmat[i,j]*100,1)), color = 'white', va = 'center', ha = 'center')

    if '--savefig' in sys.argv:
        fig3.savefig(os.path.splitext(sys.argv[1])[0]+'-gocomponent-targets.jpg', bbox_inches = 'tight', dpi = 300)


plt.show()
    





