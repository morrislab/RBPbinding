import numpy as np
import sys, os
import matplotlib.pyplot as plt
from Bio import Phylo

# Read in phylogenetic tree from timetree
tree = Phylo.read(sys.argv[1], "newick")

# Define size of figure
fig = plt.figure(figsize = (5,8.5), dpi = 300)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(left = False, labelleft = False)


# If labelfunc returns None then labels won't be shown as leafs
def labelfunc(x):
    return None

# Name of all leafs
leafs =  tree.get_terminals()

#print tree.get_nonterminals()
#print tree.find_clades({"name": "Homo_sapiens"})

# generate pairwise distance matrix between species
dmat = np.zeros((len(leafs),len(leafs)))
for l, leaf in enumerate(leafs):
    for g, geaf in enumerate(leafs):
#        print leaf, geaf
#        print tree.distance(leaf, geaf)
        dmat[l,g] = tree.distance(leaf, geaf)/2.


#tree.clade[0, 1].color = "blue"

Phylo.draw(tree, do_show=False, branch_labels = None, show_confidence = False, axes = ax, label_func = labelfunc)

amax = np.around(np.amax(dmat),1)
xticks = np.append(np.arange(amax,0, -1000), [0])
xticksmin = np.arange(amax-250,100, -250)
xticklabels = np.append(np.arange(0, amax, 1000), [amax])
ax.set_xticks(xticks.T)
ax.set_xticks(xticksmin, minor = True)
ax.tick_params(which = 'minor', bottom = True)
ax.set_xticklabels(xticklabels.T)
ax.set_xlim([-10,amax])
ax.grid(which = 'minor')
ax.set_ylabel('')
ax.set_xlabel('MyA')

if '--outname' in sys.argv:
    outname = sys.argv[2]
else:
    outname = sys.argv[1]+'.jpg'
fig.savefig(outname, dpi = 300, bbox_inches = 'tight')


