# design evaluation.py
import numpy as np
import sys, os
import matplotlib
if "--savefig" in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp
import scipy
from scipy.stats import pearsonr






rankfiles = sys.argv[1].split(',')
total = float(sys.argv[2])
names = sys.argv[3].split(',')
ranks = []
dists = []
positions = []
for r, rfile in enumerate(rankfiles):
    ranks.append(np.genfromtxt(rfile)[:,1])
    dists.append(np.genfromtxt(rfile)[:,2])
    posarray = np.zeros(int(total))
    unranks, unranknum = np.unique(ranks[-1].astype(int), return_counts = True)
    posarray[unranks] = unranknum
    positions.append(np.cumsum(posarray)/float(len(ranks[-1])))
random = np.cumsum(np.ones(int(total))/total)

colors = ['limegreen', 'dimgrey', 'purple', 'steelblue', 'firebrick']

'''
figd = plt.figure(figsize = (4,4))
axd = figd.add_subplot(111)
axd.plot(np.arange(total)+1,random, color = 'grey', alpha = 1, ls = '--', lw = 2.)
for i in range(len(ranks)):
    axd.scatter(np.arange(total)+1,positions[i], color = colors[i], alpha=  1, label = names[i])
axd.spines['top'].set_visible(False)
axd.spines['right'].set_visible(False)
axd.set_xlabel('Rank')
axd.set_ylabel('Number RBPs')
axd.set_xscale('log')
axd.set_yscale('log')
axd.set_ylim([np.amin(np.concatenate(positions)), 1.])
axd.set_yticks([0.01,0.02,0.05,0.1,0.2,0.5,1.])
axd.set_yticklabels(['1%', '2%', '5%', '10%', '20%', '50%', '100%'])
axd.set_xticks([1,2,3,5,10,20,50,100,200,500,1000])
axd.set_xticklabels([1,2,3,5,10,20,50,100,200,500,1000])
axd.set_xlim([0.9,100])
axd.legend()
'''
figd = plt.figure(figsize = (3.5,3.5), dpi = 300)
axd = figd.add_subplot(111)

for i in range(len(ranks)):
    axd.bar(np.arange(8),positions[i][[0,1,2,4,9,24,49,99]], color = colors[i], alpha=  1, label = names[i])
axd.bar(np.arange(8),random[[0,1,2,4,9,24,49,99]], color = 'k', alpha = 1, label = 'Random')
axd.spines['top'].set_visible(False)
axd.spines['right'].set_visible(False)
axd.set_xlabel('Rank sequence')
axd.set_ylabel('% tested specificities')
axd.set_ylim([np.amin(np.concatenate(positions)), 0.5])
axd.set_yticks([0.05,0.1,0.2,0.3,0.4,0.5])
axd.set_yticklabels(['5%', '10%', '20%', '30%', '40%', '50%'])
axd.set_xticks(np.arange(8))
axd.set_xticklabels([1,2,3,5,10,25,50,100])
axd.legend()


if '--savefig' in sys.argv:
    figname = sys.argv[sys.argv.index('--savefig')+1]
    figd.savefig(os.path.splitext(figname)[0]+os.path.splitext(figname)[1], dpi = 300, bbox_inches = 'tight')
else:
    plt.tight_layout()
    plt.show()






