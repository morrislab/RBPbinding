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




pset = np.genfromtxt(sys.argv[1], dtype = str)
folder = os.path.splitext(sys.argv[1])[0]
crow, srow = sys.argv[2].split(',')
pnames = np.genfromtxt(sys.argv[3], dtype = str)
pwmfile = sys.argv[4]

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

setprotnames = []
setgenenames = []
nameset = []
for pr in pset:
    setprotnames.append(pnames[list(pnames[:, 0]).index(pr), 3])
    setgenenames.append(pnames[list(pnames[:, 0]).index(pr), 2].split('.')[0])
    nameset.append(pnames[list(pnames[:, 0]).index(pr), 3]+' : \n['+pnames[list(pnames[:, 0]).index(pr), 2].split('.')[0]+']')
    
pwsort = []
for gname in pset:
    pwsort.append(pwms[list(pwmname).index(gname)])
pwms = pwsort


gofile = np.load(sys.argv[5])
ungenes = gofile['genes']

gopmat = gofile['gopmat']
gofmat = gofile['gofmat']
gocmat = gofile['gocmat']
gomat = np.concatenate([gopmat, gofmat, gocmat], axis = 1)

gotermp = gofile['gotermp']

gotermf = gofile['gotermf']

gotermc = gofile['gotermc']

ungotype = np.concatenate([['P' for i in range(len(gotermp))], ['F' for i in range(len(gotermf))], ['C' for i in range(len(gotermc))]])
ungoterm = np.concatenate([gotermp, gotermf, gotermc])

genes = []
goterm = []
gotype = []
for g, ung in enumerate(ungenes):
    has = np.where(gomat[g] > 0)[0]
    for h in has:
        genes.append(ung)
        goterm.append(ungoterm[h])
        gotype.append(ungotype[h])



'''
gofile = open('ATH_GO_GOSLIM.txt', 'r').readlines()
genes = []
goterm = []
goid = []
gotype = []
for l, line in enumerate(gofile):
    if line[0] != '!':
        line = line.strip().split('\t')
        genes.append(line[0])
        goterm.append(line[4])
        goid.append(line[5])
        gotype.append(line[7])
goid = np.array(goid)
'''
genes = np.array(genes)
goterm = np.array(goterm)
gotype = np.array(gotype)



gop = []
gof = []
goc = []
out = open(os.path.splitext(sys.argv[1])[0]+'-goterm.txt', 'w')
for p, prot in enumerate(setgenenames):

    pind = np.where(genes == prot)[0]
    gop.append(np.unique(goterm[pind[gotype[pind] == 'P']]))
    gof.append(np.unique(goterm[pind[gotype[pind] == 'F']]))
    goc.append(np.unique(goterm[pind[gotype[pind] == 'C']]))
    
    if (len(gop[-1]) > 0) or (len(gof[-1]) > 0) or (len(goc[-1]) > 0):
        print '\n', nameset[p]
    if len(gop[-1]) > 0:
        print 'P'
        print '\n'.join(gop[-1])
    if len(gof[-1]) > 0:
        print 'F'
        print '\n'.join(gof[-1])
    if len(goc[-1]) > 0:
        print 'C'
        print '\n'.join(goc[-1])
    
    out.write(' '.join(prot)+'\t"'+';'.join(gop[-1])+'"\t"'+';'.join(goc[-1])+'"\t"'+';'.join(gof[-1])+'"\n')

gopunique = np.unique(np.concatenate(gop))
gofunique = np.unique(np.concatenate(gof))
gocunique = np.unique(np.concatenate(goc))



gopmat = np.zeros((len(nameset), len(gopunique)))
gofmat = np.zeros((len(nameset), len(gofunique)))
gocmat = np.zeros((len(nameset), len(gocunique)))

for p, prot in enumerate(setgenenames):
    gopmat[p,np.isin(gopunique, gop[p])] = 1.
    gofmat[p,np.isin(gofunique, gof[p])] = 1.
    gocmat[p,np.isin(gocunique, goc[p])] = 1.


tcolmap = ListedColormap(colors)



gopmask = np.where(np.sum(gopmat, axis = 1) > 0)[0]
fig = plt.figure(figsize=(2+len(gopmat[0])*0.55,len(gopmat)*0.55), dpi = 100)

aden = fig.add_subplot(261)
aden.set_position([0.1, 0.1, 0.1 , 0.8])
aden.spines['top'].set_visible(False)
aden.spines['bottom'].set_visible(False)
aden.spines['left'].set_visible(False)
aden.spines['right'].set_visible(False)
aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
Z = linkage(cdist(gopmat[gopmask], gopmat[gopmask], 'euclidean'), 'single')
dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
sorting = dn['leaves']

pheight = 0.8/float(len(sorting))
pwidth = 2.5*pheight

adent = fig.add_subplot(262)
adent.set_position([0.2+pwidth+pheight, 0.9, 0.575 , pheight])
adent.spines['top'].set_visible(False)
adent.spines['bottom'].set_visible(False)
adent.spines['left'].set_visible(False)
adent.spines['right'].set_visible(False)
adent.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
Zg = linkage(cdist(gopmat[gopmask].T, gopmat[gopmask].T, 'euclidean'), 'single')
dng = dendrogram(Zg, orientation = 'top', color_threshold=0, above_threshold_color='k', ax = adent)
gsorting = dng['leaves']
sorting = gopmask[sorting]


for n, s in enumerate(sorting):
    print n
    axp = fig.add_subplot(len(sorting), 1, n+1)
    axp.set_position([0.2, 0.1+pheight*(n+0.1), pwidth, pheight*0.8])
    plotpwm(pwms[s], axp)
    

# imshow for rbpcolors
axt = fig.add_subplot(265)
axt.set_position([0.2+pwidth, 0.1, pheight, 0.8])
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
ax.set_position([0.2+pwidth +pheight, 0.1, 0.575, 0.8])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(which = 'both', labelleft = False, left = False, labelright = True, bottom = False)
im = ax.imshow(gopmat[sorting][:, gsorting], origin = 'lower', cmap = 'Greys', aspect = 'auto', vmin = 0, vmax = 1)
ax.set_yticks(np.arange(len(sorting)+1)-0.5, minor = True)
ax.set_xticks( np.arange(len(gopunique)+1)-0.5, minor = True)
ax.grid(color='silver', linewidth=0.5, which = 'minor')
ax.set_yticks(np.arange(len(sorting)))
ax.set_yticklabels(np.array(nameset)[sorting])

ax.set_xticks(np.arange(len(gopunique))+0.25)
ax.set_xticklabels(np.array(gopunique)[gsorting], va = 'top', ha = 'right', rotation = 60)


fig.savefig(os.path.splitext(sys.argv[1])[0]+'-goprocess.jpg', bbox_inches = 'tight', dpi = 300)






fig2 = plt.figure(figsize=(2+len(gofmat[0])*0.55,len(gofmat)*0.55), dpi = 100)

aden = fig2.add_subplot(261)
aden.set_position([0.1, 0.1, 0.1 , 0.8])
aden.spines['top'].set_visible(False)
aden.spines['bottom'].set_visible(False)
aden.spines['left'].set_visible(False)
aden.spines['right'].set_visible(False)
aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
Z = linkage(cdist(gofmat, gofmat, 'euclidean'), 'single')
dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
sorting = dn['leaves']

pheight = 0.8/float(len(nameset))
pwidth = 2.5*pheight

adent = fig2.add_subplot(262)
adent.set_position([0.2+pwidth+pheight, 0.9, 0.575 , pheight])
adent.spines['top'].set_visible(False)
adent.spines['bottom'].set_visible(False)
adent.spines['left'].set_visible(False)
adent.spines['right'].set_visible(False)
adent.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
Zg = linkage(cdist(gofmat.T, gofmat.T, 'euclidean'), 'single')
dng = dendrogram(Zg, orientation = 'top', color_threshold=0, above_threshold_color='k', ax = adent)
gsorting = dng['leaves']



for n, s in enumerate(sorting):
    print n
    axp = fig2.add_subplot(len(nameset), 1, n+1)
    axp.set_position([0.2, 0.1+pheight*(n+0.1), pwidth, pheight*0.8])
    plotpwm(pwms[s], axp)
    

# imshow for rbpcolors
axt = fig2.add_subplot(265)
axt.set_position([0.2+pwidth, 0.1, pheight, 0.8])
axt.spines['top'].set_visible(False)
axt.spines['bottom'].set_visible(False)
axt.spines['left'].set_visible(False)
axt.spines['right'].set_visible(False)
axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False, labelbottom = True, labeltop = False)
imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
axt.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
axt.grid(color='silver', linewidth=0.5, which = 'minor')
axt.set_xticks([0.1])
axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'top', ha = 'right')
for i in range(len(sorting)):
    axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center')


ax = fig2.add_subplot(263)
ax.set_position([0.2+pwidth+pheight, 0.1, 0.575, 0.8])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(which = 'both', labelleft = False, left = False, labelright = True, bottom = False)
im = ax.imshow(gofmat[sorting][:, gsorting], origin = 'lower', cmap = 'Greys', aspect = 'auto', vmin = 0, vmax = 1)
ax.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
ax.set_xticks( np.arange(len(gofunique)+1)-0.5, minor = True)
ax.grid(color='silver', linewidth=0.5, which = 'minor')
ax.set_yticks(np.arange(len(nameset)))
ax.set_yticklabels(np.array(nameset)[sorting])

ax.set_xticks(np.arange(len(gofunique))+0.25)
ax.set_xticklabels(np.array(gofunique)[gsorting], va = 'top', ha = 'right', rotation = 60)



fig2.savefig(os.path.splitext(sys.argv[1])[0]+'-gofunction.jpg', bbox_inches = 'tight', dpi = 300)










fig3 = plt.figure(figsize=(2.+len(gocmat[0])*0.55,len(gocmat)*0.55), dpi = 100)

aden = fig3.add_subplot(261)
aden.set_position([0.1, 0.1, 0.1 , 0.8])
aden.spines['top'].set_visible(False)
aden.spines['bottom'].set_visible(False)
aden.spines['left'].set_visible(False)
aden.spines['right'].set_visible(False)
aden.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
Z = linkage(cdist(gocmat, gocmat, 'euclidean'), 'single')
dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = aden)
sorting = dn['leaves']


pheight = 0.8/float(len(nameset))
pwidth = 2.5*pheight
adent = fig3.add_subplot(262)
adent.set_position([0.2+pwidth+pheight, 0.9, 0.575 , pheight])
adent.spines['top'].set_visible(False)
adent.spines['bottom'].set_visible(False)
adent.spines['left'].set_visible(False)
adent.spines['right'].set_visible(False)
adent.tick_params(bottom = False, labelbottom = False, left = False, labelleft = False)
Zg = linkage(cdist(gocmat.T, gocmat.T, 'euclidean'), 'single')
dng = dendrogram(Zg, orientation = 'top', color_threshold=0, above_threshold_color='k', ax = adent)
gsorting = dng['leaves']



for n, s in enumerate(sorting):
    print n
    axp = fig3.add_subplot(len(nameset), 1, n+1)
    axp.set_position([0.2, 0.1+pheight*(n+0.1), pwidth, pheight*0.8])
    plotpwm(pwms[s], axp)
    

# imshow for rbpcolors
axt = fig3.add_subplot(265)
axt.set_position([0.2+pwidth, 0.1, pheight, 0.8])
axt.spines['top'].set_visible(False)
axt.spines['bottom'].set_visible(False)
axt.spines['left'].set_visible(False)
axt.spines['right'].set_visible(False)
axt.tick_params(which = 'both', labelleft = False, left = False, bottom = False, labelbottom = True, labeltop = False)
imt= axt.imshow(rbpcolors[sorting].reshape(-1,1), origin = 'lower', cmap = tcolmap, aspect = 'auto', vmin = 0, vmax = 2)
axt.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
axt.grid(color='silver', linewidth=0.5, which = 'minor')
axt.set_xticks([0.1])
axt.set_xticklabels(['Stabilization'], rotation = 60, va = 'top', ha = 'right')
for i in range(len(sorting)):
    axt.text(0., i, str(np.around(float(zscores[sorting[i]]),1)), color = 'white', va = 'center', ha = 'center')


ax = fig3.add_subplot(263)
ax.set_position([0.2+pwidth+pheight, 0.1, 0.575, 0.8])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(which = 'both', labelleft = False, left = False, labelright = True, bottom = False)
im = ax.imshow(gocmat[sorting][:, gsorting], origin = 'lower', cmap = 'Greys', aspect = 'auto', vmin = 0, vmax = 1)
ax.set_yticks(np.arange(len(nameset)+1)-0.5, minor = True)
ax.set_xticks( np.arange(len(gocunique)+1)-0.5, minor = True)
ax.grid(color='silver', linewidth=0.5, which = 'minor')
ax.set_yticks(np.arange(len(nameset)))
ax.set_yticklabels(np.array(nameset)[sorting])

ax.set_xticks(np.arange(len(gocunique))+0.25)
ax.set_xticklabels(np.array(gocunique)[gsorting], va = 'top', ha = 'right', rotation = 60)


fig3.savefig(os.path.splitext(sys.argv[1])[0]+'-gocomponent.jpg', bbox_inches = 'tight', dpi = 300)


plt.show()
    





