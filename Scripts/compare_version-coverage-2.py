## coverage.py
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import glob
import pandas as pd


statsfile = sys.argv[1]
stats = np.genfromtxt(statsfile, dtype = str, delimiter = ',', skip_header = 1)
statvals = open(statsfile, 'r').readline().strip().split(', ')[1:]
statspecies = stats[:, 0][1:]
stats = stats[1:][:,1:].astype(float)

statsfile2 = sys.argv[2]
stats2 = np.genfromtxt(statsfile2, dtype = str, delimiter = ',', skip_header = 1)
statvals2 = open(statsfile2, 'r').readline().strip().split(', ')[1:]
statspecies2 = stats2[:, 0][1:]
stats2 = stats2[1:][:,1:].astype(float)


if '--savefig' in sys.argv:
    plotname = sys.argv[sys.argv.index('--savefig')+1]
    plotfmt = plotname.split('.')[-1]


if '--specieslist' in sys.argv:
    lspecies = sys.argv[sys.argv.index('--specieslist')+1]
    if ',' in lspecies:
        mspec = lspecies.split(',')
    else:
        mspec = np.genfromtxt(lspecies, dtype =str)

else:
    print 'Provide list of species'
    sys.exit()
    
    
king = False
if '--kingdoms' in sys.argv:
    king = True
    kingdoms = np.genfromtxt(sys.argv[sys.argv.index('--kingdoms')+1], dtype = str)
    if ~np.array_equal(np.sort(kingdoms[:,0]), mspec):
        resort = []
        for m, sp in enumerate(mspec):
            if sp in kingdoms[:,0]:
                resort.append([m, list(kingdoms[:,0]).index(sp)])
            else:
                print sp, 'not in kingdomfile'
        resort = np.array(resort).T
        mspec = mspec[resort[0]]
        kingdoms = kingdoms[resort[1]]


#Nmeas, Nrectot, Nrncmpt, Nmotif
inputs = [statvals.index('Nmeas'), statvals.index('Nrectot'), statvals.index('Nrncmptjple'), statvals.index('Nmotifjple')]
if '--comparefull' in sys.argv:
    inputs2 = [statvals2.index('Nmeas'), statvals2.index('Nrectot'), statvals2.index('Nrncmptjple'), statvals2.index('Nmotifjple')]
elif '--comparetoid' in sys.argv:
    inputs2 = [statvals2.index('Nmeas'), statvals2.index('Nrecid'), statvals2.index('Nrncmptid'), statvals2.index('Nmotifid')]
    
else:
    print 'define comparison'
    sys.exit()

statsort = []
for m, sp in enumerate(mspec):
    statsort.append(list(statspecies).index(sp))
stats = stats[statsort]
statspecies = statspecies[statsort]

statsort2 = []
for m, sp in enumerate(mspec):
    statsort2.append(list(statspecies2).index(sp))
stats2 = stats2[statsort2]
statspecies2 = statspecies2[statsort2]

stats = stats[:, inputs]
stats2 = stats2[:, inputs2]




if king:
    kings = np.unique(kingdoms[:,1])
    print kings

    emeans = np.mean(stats, axis = 0)
    emeans2 = np.mean(stats2, axis = 0)
    
    kmeans = np.zeros((len(kings), 4))
    kmeans2 = np.zeros((len(kings), 4))
    
    for k, king in enumerate(kings): 
        kmeans[k] = np.mean(stats[kingdoms[:,1]==king], axis = 0)
        kmeans2[k] = np.mean(stats2[kingdoms[:,1]==king], axis = 0)
    
    stats = kmeans
    stats2 = kmeans2
    mspec = kings
    if '--sortkingdom' in sys.argv:
        sortk = np.array(sys.argv[sys.argv.index('--sortkingdom')+1].split(','), dtype = int)
        stats = stats[sortk]
        stats2 = stats2[sortk]
        mspec = mspec[sortk]


if '--figsize' in sys.argv:
    xs, ys = sys.argv[sys.argv.index('--figsize')+1].split(',')
    ys = float(ys)
    xs = float(xs)
else:
    ys = 0.25*len(mspec)
    xs = 3.

for m, ms in enumerate(mspec):
    mspec[m] = ms.replace('_', ' ')

    
fig = plt.figure(figsize = (xs,ys))
ax = fig.add_subplot(111)
ax.set_xlabel('Number Motifs/Motifclusters')
ax.patch.set_facecolor('whitesmoke')
ax.barh(np.arange(1, len(mspec)+1, 1)+0.2, stats[:,2], alpha=1., color='red', label = 'Motifs\nCisBP2.0', height = 0.35)
ax.barh(np.arange(1, len(mspec)+1, 1)+0.2, stats2[:,2], alpha=1., color='darkgrey', label = 'Motifs\nCisBP0.6', height = 0.35)

ax.barh(np.arange(1, len(mspec)+1, 1)-0.2, stats[:,3], alpha=1., color='chocolate', label = 'Motif clusters\nCisBP2.0', height = 0.35)
ax.barh(np.arange(1, len(mspec)+1, 1)-0.2, stats2[:,3], alpha=1., color='tan', label = 'Motif cluster\nCisBP0.6', height = 0.35)

if king:
    ax.barh(np.array([len(mspec)+1.5])+0.2, [emeans[2]], alpha=1., color='red', height = 0.35)
    ax.barh(np.array([len(mspec)+1.5])+0.2, [emeans2[2]], alpha=1., color='darkgrey', height = 0.35)

    ax.barh(np.array([len(mspec)+1.5])-0.2, [emeans[3]], alpha=1., color='chocolate', height = 0.35)
    ax.barh(np.array([len(mspec)+1.5])-0.2, [emeans2[3]], alpha=1., color='tan',  height = 0.35)

    ax.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
    ax.set_yticklabels(np.append(mspec, ['Eukaryotes']))

else:
    ax.set_yticks(np.arange(1., len(mspec)+1., 1))
    ax.set_yticklabels(mspec)
minor_xticks = np.arange(0,np.amax(stats[:, 2]), 10, dtype = int)
#minor_xticks = np.linspace(0.05,1.05,11)
ax.set_xticks(minor_xticks)
ax.set_xticklabels(minor_xticks, rotation = 60)
ax.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
#ax.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_linewidth(2.)
if '--legend' in sys.argv:
    ax.legend()


fig2 = plt.figure(figsize = (xs,ys))
ax2 = fig2.add_subplot(111)
ax2.set_xlabel('Number RBPs')
ax2.patch.set_facecolor('whitesmoke')
ax2.barh(np.arange(1, len(mspec)+1, 1), stats[:,1], alpha=1., color='red', label = 'Reconstructed\nRBPs CisBP2.0')
ax2.barh(np.arange(1, len(mspec)+1, 1), stats2[:,1], alpha=1., color='dimgrey', label = 'Reconstruced\nRBPs CisBP0.6')
ax2.barh(np.arange(1, len(mspec)+1, 1), stats[:,0], alpha=1., color='cornflowerblue', label = 'Measured\nRBPs CisBP2.0')
ax2.barh(np.arange(1, len(mspec)+1, 1), stats2[:,0], alpha=1., color='black', label = 'Measured\nRBPs CisBP0.6')


if king:
    # eukaryotes
    ax2.barh(np.array([len(mspec)+1.5]), [emeans[1]], alpha=1., color='red')
    ax2.barh(np.array([len(mspec)+1.5]), [emeans2[1]], alpha=1., color='dimgrey')
    ax2.barh(np.array([len(mspec)+1.5]), [emeans[0]], alpha=1., color='cornflowerblue')
    ax2.barh(np.array([len(mspec)+1.5]), [emeans2[0]], alpha=1., color='black')
    
    ax2.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
    ax2.set_yticklabels(np.append(mspec, ['Eukaryotes']))

else:
    ax2.set_yticks(np.arange(1., len(mspec)+1., 1))
    ax2.set_yticklabels(mspec)

major_xticks = np.arange(0, np.amax(stats[:,1]), 20, dtype = int)
#minor_xticks = np.linspace(0.05,1.05,11)
ax2.set_xticks(major_xticks)
ax2.set_xticklabels(major_xticks, rotation =60)
ax2.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
#ax.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
#ax.spines['bottom'].set_linewidth(2.)
if '--legend' in sys.argv:
    ax2.legend()


fig3 = plt.figure(figsize = (xs,ys))
ax3 = fig3.add_subplot(111)
ax3.set_xlabel('Number Motifs/Motifclusters')
ax3.patch.set_facecolor('whitesmoke')
ax3.barh(np.arange(1, len(mspec)+1, 1), stats[:,2], alpha=1., color='red', label = 'Motifs CisBP2.0')
ax3.barh(np.arange(1, len(mspec)+1, 1), stats2[:,2], alpha=1., color='darkgrey', label = 'Motifs CisBP0.6')

if king:
    ax3.barh(np.array([len(mspec)+1.5]), [emeans[2]], alpha=1., color='red')
    ax3.barh(np.array([len(mspec)+1.5]), [emeans2[2]], alpha=1., color='darkgrey')

    ax3.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
    ax3.set_yticklabels(np.append(mspec, ['Eukaryotes']))

else:
    ax3.set_yticks(np.arange(1., len(mspec)+1., 1))
    ax3.set_yticklabels(mspec)
minor_xticks = np.arange(0,np.amax(stats[:, 2]), 10, dtype = int)
#minor_xticks = np.linspace(0.05,1.05,11)
ax3.set_xticks(minor_xticks)
ax3.set_xticklabels(minor_xticks, rotation = 60)
ax3.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
#ax.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['left'].set_visible(False)
#ax.spines['bottom'].set_linewidth(2.)
if '--legend' in sys.argv:
    ax3.legend()

fig4 = plt.figure(figsize = (xs,ys))
ax4= fig4.add_subplot(111)
ax4.set_xlabel('Number Motifs/Motifclusters')
ax4.patch.set_facecolor('whitesmoke')
ax4.barh(np.arange(1, len(mspec)+1, 1), stats[:,3], alpha=1.,  color='chocolate', label = 'Motif clusters CisBP2.0')
ax4.barh(np.arange(1, len(mspec)+1, 1), stats2[:,3], alpha=1., color='tan', label = 'Motif clusters CisBP0.6')

if king:
    ax4.barh(np.array([len(mspec)+1.5]), [emeans[3]], alpha=1.,  color='chocolate')
    ax4.barh(np.array([len(mspec)+1.5]), [emeans2[3]], alpha=1., color='tan')

    ax4.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
    ax4.set_yticklabels(np.append(mspec, ['Eukaryotes']))

else:
    ax4.set_yticks(np.arange(1., len(mspec)+1., 1))
    ax4.set_yticklabels(mspec)
minor_xticks = np.arange(0,np.amax(stats[:, 3]), 10, dtype = int)
#minor_xticks = np.linspace(0.05,1.05,11)
ax4.set_xticks(minor_xticks)
ax4.set_xticklabels(minor_xticks, rotation = 60)
ax4.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
#ax.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['left'].set_visible(False)
#ax.spines['bottom'].set_linewidth(2.)
if '--legend' in sys.argv:
    ax4.legend()


if '--savefig' in sys.argv:
    fig.savefig(os.path.splitext(plotname)[0]+'_motifreconst.'+plotfmt, dpi = 300, bbox_inches = 'tight')
    fig2.savefig(os.path.splitext(plotname)[0]+'_rbpreconst.'+plotfmt, dpi = 300, bbox_inches = 'tight')
    fig3.savefig(os.path.splitext(plotname)[0]+'_motifonlyreconst.'+plotfmt, dpi = 300, bbox_inches = 'tight')
    fig4.savefig(os.path.splitext(plotname)[0]+'_clusterreconst.'+plotfmt, dpi = 300, bbox_inches = 'tight')
else:
    fig4.tight_layout()
    fig3.tight_layout()
    fig2.tight_layout()
    fig.tight_layout()
    plt.show() 




