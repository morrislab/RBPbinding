import numpy as np
import sys, os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.cluster import hierarchy
import logomaker as lm

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
for p,pr in enumerate(pset):
    if '--double' in sys.argv:
        setprotname = pnames[list(pnames[:, 0]).index(pr.split('--')[0]), 3]
        setgenename = pnames[list(pnames[:, 0]).index(pr.split('--')[1]), 3].split('.')[0]
        nameset.append(setprotname+': '+zscores[p][:4]+'\n-'+setgenename)
    else:
        setprotname = pnames[list(pnames[:, 0]).index(pr), 3]
        setgenename = pnames[list(pnames[:, 0]).index(pr), 2].split('.')[0]
        nameset.append(setprotname+': '+zscores[p][:4]+'\n['+setgenename+']')

bexpsort = []
pwsort = []
for gname in pset:
    bexpsort.append(list(rbpnames).index(gname))
    if '--double' in sys.argv:
        pwsort.append(np.concatenate([pwms[list(pwmname).index(gname.split('--')[0])], np.ones((2,4))*0.25, pwms[list(pwmname).index(gname.split('--')[1])]],axis =0))
    else:
        pwsort.append(pwms[list(pwmname).index(gname)])
pwms = pwsort
correlationmat = correlationmat[bexpsort]
bindingmat = bindingmat[bexpsort]



bind = [1]
sigcor = 0.2


genecorrelation = []
genecornum = []
genecornumcor = []
ngencontrol = []


numberset = []
for p, correlation in enumerate(correlationmat):
    ngencontrol.append(correlation[~bindingmat[p]])
    genecorrelation.append(correlation[bindingmat[p]])
    genecornum.append(np.sum(bindingmat[p]))
    genecornumcor.append(np.sum(np.sign(float(zscores[p]))*correlation[bindingmat[p]] > sigcor))
    
    numberset.append(str(int(genecornum[-1]))+'/'+str(int(genecornumcor[-1])))
    

fig = plt.figure(figsize=(8,len(genecorrelation)*0.7))
ax = fig.add_subplot(121)
ax.set_position([0.1, 0.1, 0.5, 0.8])

d2 = ax.boxplot(list(genecorrelation), positions = np.arange(len(pset))+0.25, widths = 0.25, patch_artist = True,vert = False, showfliers = False, boxprops=dict(facecolor='slateblue', color='k'))
d3 = ax.boxplot(list(ngencontrol),positions = np.arange(len(pset)), widths = 0.25, patch_artist = True, vert= False, showfliers = False, boxprops=dict(facecolor='silver', color='k'))



ax.set_yticks(np.arange(len(nameset))+0.25)
ax.set_yticklabels(nameset, rotation = 0)
ax.set_ylim([-.25, len(nameset)-0.25])
ax.plot([0.,0.], [0.,len(nameset)-0.25], ls = '-', c = 'grey', alpha = 0.4)
ax.plot([-sigcor, -sigcor], [0.,len(nameset)-0.25], ls = '--', c = 'grey', alpha = 0.4)
ax.plot([sigcor, sigcor], [0.,len(nameset)-0.25], ls = '--', c = 'grey', alpha = 0.4)
ax.set_xlabel('Pearson R')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim([-1.05,1.05])

ylim = ax.get_ylim()

pheight = 0.8/float(len(nameset))
pwmwidth = 0.1 *float(len(pwms[0]))/9.
print pwmwidth

axbar = fig.add_subplot(122)
axbar.set_position([0.6+pwmwidth + 0.01,0.1,0.2,0.8])

axbar.barh(np.arange(len(pset)), genecornum, height = 0.25, color = 'slateblue', edgecolor = 'grey')
axbar.barh(np.arange(len(pset)), genecornumcor, height = 0.25, color = 'lightslategrey', edgecolor = 'grey')
axbar.spines['right'].set_visible(False)
axbar.spines['top'].set_visible(False)
axbar.spines['left'].set_visible(False)
axbar.set_xlabel('Number\ntargets')
axbar.tick_params(left = False, labelleft = False, labelright = True)
axbar.set_yticks(np.arange(len(pset))+0.25)
axbar.set_yticklabels(numberset)
axbar.set_ylim(ylim)


for n in range(len(nameset)):
    axp = fig.add_subplot(len(nameset), 1, n+1)
    axp.set_position([0.6, 0.1+pheight*n+pheight*0.1, pwmwidth, pheight*0.8])
    plotpwm(pwms[n], axp)

if '--savefig' in sys.argv:
    fig.savefig(outname+'_correlation_distribution_bound_unbound.jpg', dpi = 300, bbox_inches = 'tight')
else:
    plt.show()



