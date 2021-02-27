import numpy as np
import sys, os
import matplotlib
if "--savefig" in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
from scipy.stats import wilcoxon
import matplotlib.cm as cm
import seaborn as sb
from functools import reduce 
from scipy.stats import gaussian_kde


if '--savefig' in sys.argv:
    plotname = sys.argv[sys.argv.index('--savefig')+1]
    plotfmt = plotname.rsplit('.')[-1]
    


if "--files" in sys.argv:
    files = []
    colums = []
    for i in range(2):
        if "--definecolums" in sys.argv:
            files.append(sys.argv[sys.argv.index("--files")+1+i])
            colums.append(int(sys.argv[sys.argv.index("--definecolums")+1+i]))
        elif "--defineallcolums" in sys.argv:
            files.append(sys.argv[sys.argv.index("--files")+1+i])
            colums.append(int(sys.argv[sys.argv.index("--defineallcolums")+1]))
        else:
            files.append(sys.argv[sys.argv.index("--files")+1+i])
            colums.append(-2)
if '--namerow' in sys.argv:
    nameline = int(sys.argv[sys.argv.index('--namerow')+1])
else:
    nameline = 0
    
correlations = []
sets = []
protnames = []
for i, infile in enumerate(files):
    #print infile
    protnames.append(np.genfromtxt(infile, dtype = str)[:,nameline])
    #print colums[i]
    corr = np.nan_to_num(np.genfromtxt(infile)[:,colums[i]])
    #print corr
    correlations.append(corr)
    sets.append(infile)
#print np.shape(protnames)
protnames = np.array(protnames)
names = protnames[0]
correlations = np.array(correlations)
set = np.array(sets)



if '--minlength' in sys.argv:
    minlen = int(sys.argv[sys.argv.index('--minlength')+1])
    #print minlen
    keep = []
   # print correlations[-1]
    for c, corr in enumerate(correlations):
        #print sets[c], len(np.unique(protnames[c]))
        if len(np.unique(protnames[c])) >= minlen:
            keep.append(c)
        #else:
            #print sets[c], len(np.unique(protnames[c])), 'excluded'
    protnames = protnames[keep]
    #print np.shape(protnames)
    correlations = correlations[keep]
    if len( keep) < len(sets):
        print 'ONLY', len(keep), 'kept'
    sets = np.array(sets)[keep]



intersec = reduce(np.intersect1d, (protnames[:]))
print 'Common set of proteins', len(intersec)
fincorrelation = np.zeros((len(correlations),len(intersec)))
for j, prot in enumerate(intersec):
    for i in range(len(protnames)):
        fincorrelation[i,j] = correlations[i][np.where(protnames[i] == prot)[0][0]] 
names = intersec
correlations = fincorrelation[:]
newpnames = []
for i in range(len(protnames)):
    newpnames.append( names)
protnames = np.array(newpnames)[0]


if "--setnames" in sys.argv:
    # assing shorter names to the files
    for i in range(len(files)):
        sets[i]=sys.argv[sys.argv.index("--setnames")+1+i]


if '--proteinset' in sys.argv:
    pnlist = np.genfromtxt(sys.argv[sys.argv.index('--proteinset')+1], dtype = str)[:, 0]
    fincorrelation = []
    nnames = []
    print 'Proteinlist given'
    print 'Before', len(correlations[0])
    #print pnlist
    for i in range(len(protnames)):
        #print protnames[i]
        ncorrel = []
        keepind = []
        for j, prot in enumerate(pnlist):
            if prot in protnames[i]:
                ki = list(protnames[i]).index(prot)
                keepind.append(ki)
                ncorrel.append(correlations[i][ki])
            #else:
                #print prot
        nnames.append(np.array(protnames[i])[keepind])
        fincorrelation.append(ncorrel)
    protnames = nnames
    names = nnames[0]
    correlations = np.array(fincorrelation)
    plotnames = plotnames[keepind]
    print 'After', len(correlations[0])
    


if '--rescale' in sys.argv:
    scalingfile = sys.argv[sys.argv.index('--rescale')+1]
    scaleline = sys.argv[sys.argv.index('--rescale')+2]

    sfile = np.genfromtxt(scalingfile, dtype = str)

    snames = sfile[:, 0]
    scaling = sfile[:, int(scaleline)].astype(float)

    
    for c in range(len(correlations)):
        rescaled = []
        for p in range(len(protnames)):
            if protnames[p] in snames:
                rescaled.append(max(-1,min(1., correlations[c, p]/scaling[list(snames).index(protnames[p])])))
        correlations[c] = rescaled
            
    
if '--filterprots' in sys.argv:
    filterfile = sys.argv[sys.argv.index('--filterprots')+1]
    fline = sys.argv[sys.argv.index('--filterprots')+2]
    fcut = float(sys.argv[sys.argv.index('--filterprots')+3])
    fdir = sys.argv[sys.argv.index('--filterprots')+4]
    
    ffile = np.genfromtxt(filterfile, dtype = str)

    fnames = ffile[:, 0]
    fcaling = ffile[:, int(fline)].astype(float)

    keepprot = []
    for p in range(len(protnames)):
        if protnames[p] in fnames:
            keepprot.append(fcaling[list(fnames).index(protnames[p])])
    if fdir == 'greater':
        keepmask = np.array(keepprot) >= fcut
    else:
        keepmask = np.array(keepprot) <= fcut
    protnames = protnames[keepmask]
    correlations = np.array(correlations)[:, keepmask]
    



if '--plotlim' in sys.argv:
    threshold = sys.argv[sys.argv.index('--plotlim')+1]
    if threshold[-1].isdigit() or threshold[-1] == '.':
        threshold = np.array(threshold.split(','), dtype = float)
        lim = [threshold[0], threshold[-1]]
    else:
        lim = [np.amin(correlations)-0.025*np.amin(correlations), np.amax(correlations)+0.025*np.amax(correlations)]

tresh = False
if '--threshold' in sys.argv:
    tresh = True
    threshold = float(sys.argv[sys.argv.index('--threshold')+1])
    
if '--colorfile' in sys.argv:
    colfile = sys.argv[sys.argv.index('--colorfile')+1]
    creverse = False
    if colfile[0] == '-':
        creverse = True
        colfile = colfile[1:]
    if ':' in colfile:
        colcol = int(colfile.rsplit(':',1)[-1])
        colfile = colfile.rsplit(':',1)[0]
    else:
        colcol = -1
    colfile = np.genfromtxt(colfile, dtype = str)
    cpnames = colfile[:,0]
    cps = colfile[:,colcol].astype(float)
    if creverse:
        cps = 1.-cps
    z = []
    for name in protnames:
        z.append(cps[list(cpnames).index(name)])
    z = np.array(z)
elif '--colordensity' in sys.argv:
    xy = np.vstack([correlations[0],correlations[1]])
    z = np.log(1.+gaussian_kde(xy)(xy))
else:
    z = np.ones(len(correlations[0]))


sze = np.ones(len(correlations[0]))*60
if '--sizefile' in sys.argv:
    sizefile = sys.argv[sys.argv.index('--sizefile')+1]
    creverse = False
    if sizefile[0] == '-':
        creverse = True
        sizefile = sizefile[1:]
    if ':' in sizefile:
        sizecol = int(sizefile.rsplit(':',1)[-1])
        sizefile = sizefile.rsplit(':',1)[0]
    else:
        sizecol = -1
    sizefile = np.genfromtxt(sizefile, dtype = str)
    
    cpnames = sizefile[:,0]
    cps = sizefile[:,sizecol].astype(float)
    if creverse:
        cps = 1.-cps
        #print cps
    sze = []
    for name in protnames:
        sze.append(cps[list(cpnames).index(name)])
    sze = np.array(sze)



fig = plt.figure(figsize = (3,3))
ax = fig.add_subplot(111)
ax.set_position([0.1, 0.1, 0.8, 0.8 ])
if '--colorfile' in sys.argv and '--colorbar' in sys.argv:
    cmin = float(sys.argv[sys.argv.index('--colorbar')+1])
    cmax = float(sys.argv[sys.argv.index('--colorbar')+2])
    cbar = fig.add_subplot(3,2,1)
    cbar.set_position([0.1, .95, 0.15, 0.05 ])
    cbar.imshow(np.linspace(cmin, cmax, 100).reshape(1,-1), cmap = cm.Blues_r, aspect = 'auto', vmin = cmin, vmax = 1.2*cmax)
    cbar.spines['top'].set_visible(False)
    cbar.spines['bottom'].set_visible(False)
    cbar.spines['right'].set_visible(False)
    cbar.spines['left'].set_visible(False)
    cbar.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False, labeltop = True, top = True)
    cbar.set_xticks([0, 50, 99])
    cbar.set_xticklabels(np.array([cmin, np.around((cmax-cmin)/2., 1),cmax], dtype = int), rotation = 60)
else:
    cmin = np.amin(z)
    cmax = np.amax(z)

if '--sizefile' in sys.argv and '--colorbar' in sys.argv:
    
    sbar = fig.add_subplot(3,2,2)
    cbar.set_position([0.75, .95, 0.15, 0.05 ])
    cbar.scatter([0, 0.5, 1.], np.around([np.amin(sze), (np.amin(sze)+np.amax(sze))/2.,np.amax(sze)],3), c = 'lightgrey', edgecolor = 'k', size = [np.amin(sze), (np.amin(sze)+np.amax(sze))/2.,np.amax(sze)])
    cbar.spines['top'].set_visible(False)
    cbar.spines['bottom'].set_visible(False)
    cbar.spines['right'].set_visible(False)
    cbar.spines['left'].set_visible(False)
    cbar.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False, labeltop = True, top = True)
    cbar.set_xticks([0, 0.5, 1.])
    cbar.set_xticklabels(np.around([np.amin(sze), (np.amin(sze)+np.amax(sze))/2.,np.amax(sze)],3), rotation = 60)    

label = ''
if '--label' in sys.argv:
    label = sys.argv[sys.argv.index('--label')+1]


csort = np.argsort(-z)
afi = ax.scatter(correlations[0][csort],correlations[1][csort], cmap=cm.Blues_r, c = z[csort], s = sze[csort], alpha = 0.7, vmin = cmin, vmax = cmax+0.2*cmax)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

    
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.plot([lim[0],lim[-1]], [lim[0],lim[-1]], c= 'gray')
if tresh:
    ax.plot([lim[0], lim[-1]], [threshold,threshold], c= 'gray', linestyle = '--', alpha = 0.5)
    ax.plot([threshold,threshold], [lim[0], lim[-1]], c= 'gray', linestyle = '--', alpha = 0.5)

ax.set_xlabel(label+' '+sets[0]+' (Mean: '+str(np.around(np.mean(correlations[0]),2))+')')
ax.set_ylabel(label+' '+sets[1]+' (Mean: '+str(np.around(np.mean(correlations[1]),2))+')')

#fig.tight_layout()
if '--savefig' in sys.argv:
    fig.savefig(plotname, format=plotfmt, bbox_inches = 'tight', dpi = 300)
else:
    plt.show()             
