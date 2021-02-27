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






    
rbpfile = np.load(sys.argv[1])
rbpexpression = rbpfile['pkm']
tissues = rbpfile['tissues']
rbps = rbpfile['genes']

expressionfile = np.load(sys.argv[2])
expression = expressionfile['pkm']
etissues = expressionfile['tissues']
genes = expressionfile['genes']

if not np.array_equal(tissues, etissues):
    print 'Tissues are not equal', np.array_equal(tissues, etissues)
    sys.exit()

tfile = open('SraRunTable.txt', 'r').readlines()
tissueid = []
for l, line in enumerate(tfile):
    if l> 0:
        line = line.split('\t')
        tissueid.append([line[8], line[10]])
tissueid = np.array(tissueid)
uniquetissue, tn = np.unique(tissueid[:, 1], return_index = True)
nrbpexpression = np.zeros((len(rbps), len(uniquetissue)))
nexpression = np.zeros((len(genes), len(uniquetissue)))
for u, utiss in enumerate(uniquetissue):
    tissueids = tissueid[tissueid[:,1] == utiss, 0]
    loc = np.where(np.isin(tissues, tissueids))[0]
    nrbpexpression[:, u] = np.mean(rbpexpression[:, loc], axis = 1)
    nexpression[:, u] = np.mean(expression[:, loc], axis = 1)
rbpexpression = nrbpexpression
expression = nexpression
tissues = uniquetissue
print 'mean tissue', np.shape(rbpexpression)


pset = np.genfromtxt(sys.argv[3], dtype = str)
folder = os.path.splitext(sys.argv[3])[0]
crow, srow = sys.argv[4].split(',')
pnames = np.genfromtxt(sys.argv[5], dtype = str)
pwmfile = sys.argv[6]
pwmname, pwms = readinpwm(pwmfile)
pset = pset[pset[:, int(srow)] == 'True']

colors = np.array(['grey', 'coral', 'cornflowerblue'])
rbpcolors = np.zeros(len(pset), dtype = int)
rbpcolors[pset[:, int(crow)].astype(float) > 0] = 1
rbpcolors[pset[:, int(crow)].astype(float) < 0] = 2
zscores = pset[:, int(crow)]
pset = pset[:, 0]

zsort = np.argsort(zscores.astype(float))

pset = pset[zsort]
zscores = zscores[zsort]
rbpcolors = rbpcolors[zsort]





if '--double' in sys.argv:
    rbpexpression = (rbpexpression.T-np.amin(rbpexpression, axis = 1)).T
    rbpexpression = (rbpexpression.T/np.amax(rbpexpression, axis = 1)).T
    nrbpexpression = []
    npwms = []
    for p, pr in enumerate(pset):
        prs = pr.split('--')
        nrbpexpression.append(rbpexpression[list(rbps).index(prs[0])]*rbpexpression[list(rbps).index(prs[1])])
        npwms.append(np.concatenate([pwms[list(pwmname).index(prs[0])], np.ones((2,4))*0.25, pwms[list(pwmname).index(prs[1])]],axis =0))
    rbpexpression = np.array(nrbpexpression)
    pwms = npwms
    rbps = pset
else:
    rexpsort = []
    pwsort = []
    for gname in pset:
        rexpsort.append(list(rbps).index(gname))
        pwsort.append(pwms[list(pwmname).index(gname)])
    pwms = pwsort
    rbpexpression = rbpexpression[rexpsort]
    rbps = np.array(rbps)[rexpsort]



selectivity0 = (np.amax(rbpexpression, axis = 1)-np.mean(rbpexpression, axis = 1))/np.std(rbpexpression, axis = 1)
selectivity = skew(rbpexpression, axis = 1)
exp5,exp50,exp95 = np.percentile(rbpexpression, [5, 50, 95], axis = 1)
selectivity3 = (exp95-exp50)/(exp50-exp5)


add = ''
if '--filter_selectivity' in sys.argv:
    selectmin = float(sys.argv[sys.argv.index('--filter_selectivity')+1])
    direction = sys.argv[sys.argv.index('--filter_selectivity')+2]
    if direction == 'gt':
        rbpmask = selectivity > selectmin
    else:
        rbpmask = selectivity <= selectmin
    add+= '_filts-'+direction+str(selectmin)
    rbpexpression = rbpexpression[rbpmask]
    print np.shape(rbpexpression)
    selectivity = selectivity[rbpmask]
    pset = pset[rbpmask]
    zscores = zscores[rbpmask]
    rbpcolors = rbpcolors[rbpmask]
    pwms = np.array(pwms)[rbpmask]
    rbps = np.array(rbps)[rbpmask]


if '--clustertissue' in sys.argv:
    clusterfile = sys.argv[sys.argv.index('--clustertissue')+1]
    add += '_tissucluster' + os.path.splitext(os.path.split(clusterfile)[1])[0]
    clusters = np.genfromtxt(clusterfile, dtype = str, delimiter = '\t')
    uniqueclusters, cind = np.unique(clusters[:,1], return_index = True)
    uniqueclusters = uniqueclusters[np.argsort(cind)]
    nrbpexpression = np.zeros((len(pset), len(uniqueclusters)))
    nexpression = np.zeros((len(genes), len(uniqueclusters)))
    for u, utiss in enumerate(uniqueclusters):
        loc = np.where(np.isin(tissues, clusters[clusters[:,1]==utiss,0]))[0]
        nrbpexpression[:, u] = np.mean(rbpexpression[:, loc], axis = 1)
        nexpression[:, u] = np.mean(expression[:, loc], axis = 1)
    rbpexpression = nrbpexpression
    expression = nexpression
    tissues = uniqueclusters
    print 'mean tissue', np.shape(rbpexpression)

#if '--double' in sys.argv:
    #normrbpexpression = rbpexpression
normrbpexpression = (rbpexpression.T-np.amin(rbpexpression, axis = 1)).T
normrbpexpression = (normrbpexpression.T/np.amax(normrbpexpression, axis = 1)).T

if '--tissuefilter' in sys.argv:
    clusterfile = sys.argv[sys.argv.index('--tissuefilter')+1]
    add += '_tissues' + os.path.splitext(os.path.split(clusterfile)[1])[0]   
    clusters = np.genfromtxt(clusterfile, dtype = str, delimiter = '\n')
    keep = []
    for u, utiss in enumerate(clusters):
        keep.append(list(tissues).index(utiss))
    rbpexpression = rbpexpression[:, keep]
    normrbpexpression = normrbpexpression[:, keep]
    expression = expression[:, keep]
    tissues = clusters
    if '--filter_lowrbps' in sys.argv:
        rbpmask =  np.amax(normrbpexpression, axis =1) >= 0.8
        rbpexpression = rbpexpression[rbpmask]
        normrbpexpression = normrbpexpression[rbpmask]
        print np.shape(rbpexpression)
        selectivity = selectivity[rbpmask]
        pset = pset[rbpmask]
        zscores = zscores[rbpmask]
        rbpcolors = rbpcolors[rbpmask]
        pwms = np.array(pwms)[rbpmask]
        rbps = np.array(rbps)[rbpmask]
        
selectivity0 = (np.amax(rbpexpression, axis = 1)-np.mean(rbpexpression, axis = 1))/np.std(rbpexpression, axis = 1)
selectivity = skew(rbpexpression, axis = 1)
exp5,exp50,exp95 = np.percentile(rbpexpression, [5, 50, 95], axis = 1)
selectivity3 = (exp95-exp50)/(exp50-exp5)

nameset = []
for pr in pset:
    if '--double' in sys.argv:
        setprotname = pnames[list(pnames[:, 0]).index(pr.split('--')[0]), 3]+'-'+pnames[list(pnames[:, 0]).index(pr.split('--')[1]), 3]
    else:
        setprotname = pnames[list(pnames[:, 0]).index(pr), 3]+'\\'+pnames[list(pnames[:, 0]).index(pr), 2]
    nameset.append(setprotname)
nameset = np.array(nameset)


cmapgs = np.concatenate([[[1.,1.,1.,1.] for i in range(40)], [cm.Greens(0.55) for i in range(30)], [(0.,0.1,0.03, 1.) for i in range(30)]], axis = 0)
#cmapgs = np.concatenate([[[1.,1.,1.,1.] for i in range(20)] , [cm.Greens(0.2) for i in range(20)], [cm.Greens(0.55) for i in range(20)], [cm.Greens(0.95) for i in range(20)], [(0.,0.05,0., 1.) for i in range(20)]], axis = 0)

colmap = ListedColormap(cmapgs)
tcolmap = ListedColormap(colors)

selctcolor = np.concatenate([[cm.Greys((50.-i)/100.) for i in range(50)], [[1.,1.,1.,1.] for i in range(50)], [cm.Oranges(i/50.) for i in range(50)]], axis = 0)
selctcolor = ListedColormap(selctcolor)


fig = plt.figure(figsize=(3+0.3* len(tissues),1+len(rbpexpression)*0.3), dpi = 50)
wunit = 0.025


ac = fig.add_subplot(266)
ac.set_position([0.9 -2*wunit, 0.95, 1.5*wunit ,0.05])
ac.imshow([np.linspace(0,1, 100)], origin = 'lower', cmap = colmap, aspect = 'auto')
ac.set_xticks([0,99])

ac.set_xticklabels(['Min','Max'])
#ac.set_xticklabels([0,0.25,0.5,0.75,1.])
#ac.spines['top'].set_visible(False)
#ac.spines['left'].set_visible(False)
#ac.spines['right'].set_visible(False)
#ac.spines['bottom'].set_visible(False)        
ac.tick_params(left = False, labelleft = False, right = False, labelright = False)
#ac.set_title('Top 100 overlap')
ac.set_title('RPKM')

if not '--sortedtissue' in sys.argv:
    add += '-sortedtis'
    adent = fig.add_subplot(261)
    adent.set_position([0.1+2*wunit, 0.9, 0.8 - 5*wunit , 0.1])
    if '--fullsimilarity' in sys.argv:
        add += '-fullsim'
        tissueorder = cdist(expression.T, expression.T, 'correlation')
    else:
        tissueorder = cdist(normrbpexpression.T, normrbpexpression.T, 'correlation')
    Ztissue = linkage(tissueorder, 'single')
    dnt = dendrogram(Ztissue, orientation = 'top', color_threshold=0, above_threshold_color='k', ax = adent)
    tsorting = dnt['leaves']
    adent.spines['top'].set_visible(False)
    adent.spines['bottom'].set_visible(False)
    adent.spines['left'].set_visible(False)
    adent.spines['right'].set_visible(False)
    adent.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
else:
    tsorting = np.arange(len(tissues))

aden = fig.add_subplot(262)
aden.set_position([0.1-2*wunit, 0.1, 4*wunit , 0.8])
aden.spines['top'].set_visible(False)
aden.spines['bottom'].set_visible(False)
aden.spines['left'].set_visible(False)
aden.spines['right'].set_visible(False)
aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
idmat = cdist(normrbpexpression, normrbpexpression, 'euclidean')
Z = linkage(idmat, 'average')
dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
sorting = dn['leaves']

ax = fig.add_subplot(263)
ax.set_position([0.1+2*wunit, 0.1, 0.8-5*wunit, 0.8])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(which = 'both', labelleft = False, left = False, bottom = False)
normrbpexpression = normrbpexpression[sorting][:, tsorting]
im = ax.imshow(normrbpexpression, origin = 'lower', cmap = colmap, aspect = 'auto', vmin = 0, vmax = 1)
ax.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
ax.set_xticks( np.arange(len(tissues)+1)-0.5, minor = True)
ax.grid(color='silver', linewidth=0.5, which = 'minor')
ax.set_xticks(np.arange(len(tissues))+0.25)
ax.set_xticklabels(tissues[tsorting], rotation = 60, va = 'top', ha = 'right') 

#for i in tsorting:
    #print tissues[i]

pheight = 0.8/float(len(nameset))
for n, s in enumerate(sorting):
    print n
    axp = fig.add_subplot(len(nameset), 1, n+1)
    axp.set_position([0.9-0.7*wunit, 0.1+pheight*(n+0.1), wunit*2.55555, pheight*0.8])
    plotpwm(pwms[s], axp)
    axp.tick_params(labelright = True)
    axp.set_yticks([.5])
    axp.set_yticklabels([nameset[s]])
    




axs = fig.add_subplot(264)
axs.set_position([0.9-2.8*wunit, 0.1, wunit, 0.8])
axs.spines['top'].set_visible(False)
axs.spines['bottom'].set_visible(False)
axs.spines['left'].set_visible(False)
axs.spines['right'].set_visible(False)
axs.tick_params(which = 'both', labelleft = False, left = False, bottom = False)
ims = axs.imshow(selectivity[sorting].reshape(-1,1), origin = 'lower', cmap = selctcolor, aspect = 'auto', vmin = -1.5, vmax = 1.5)
axs.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
axs.grid(color='silver', linewidth=0.5, which = 'minor')
axs.set_xticks([0.25])
axs.set_xticklabels(['Selectivity'], rotation = 60, va = 'top', ha = 'right')
for i in range(len(selectivity)):
    axs.text(0., i, str(np.around(selectivity[sorting[i]],1)), color = 'grey', va = 'center', ha = 'center', fontsize = 9)

# imshow for rbpcolors
axt = fig.add_subplot(265)
axt.set_position([0.9-1.8*wunit, 0.1, wunit, 0.8])
axt.spines['top'].set_visible(False)
axt.spines['bottom'].set_visible(False)
axt.spines['left'].set_visible(False)
axt.spines['right'].set_visible(False)
axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False)
imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
axt.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
axt.grid(color='silver', linewidth=0.5, which = 'minor')
axt.set_xticks([0.25])
axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'top', ha = 'right')
for i in range(len(selectivity)):
    axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center', fontsize = 9)

figexp = plt.figure(figsize = (5, len(nameset)*0.5))
axe = figexp.add_subplot(111)
axe.spines['top'].set_visible(False)
axe.spines['left'].set_visible(False)
axe.spines['right'].set_visible(False)
axe.set_xscale('log')
axe.boxplot(list(rbpexpression[sorting]), widths = 0.75, vert = False) #, patch_artist = True, vert = False, showfliers = True, boxprops=dict(facecolor='darkslateblue', color='k'))
axe.plot([0,0], [0, len(rbpexpression)], color = 'grey')
axe.set_yticks(np.arange(len(nameset))+1)
axe.set_yticklabels(nameset[sorting], rotation = 0)
axe.set_xlabel('RPKM')


fig.savefig(folder+add+'-ex.jpg', bbox_inches = 'tight', dpi = 300)
figexp.savefig(folder+add+'-exbbox.jpg', bbox_inches ='tight', dpi =300)
plt.show()



