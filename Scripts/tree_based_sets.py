import numpy as np
import pandas as pd
from ete3 import ClusterTree
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
import logging
import sys, os
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib.colors import ListedColormap
import seaborn as sns
from scipy.spatial.distance import cdist
# make linkage with 3 worm species
import matplotlib.cm as cm
from skbio import DistanceMatrix
from skbio.tree import nj
from Bio import Phylo
from io import StringIO


# Read in different files
def readinclustermat(clustermat, splitnames = None, namepart = 1, maxid = 100.):
    filetype = os.path.splitext(clustermat)[-1]
    assert filetype in ['.txt', '.dat', '.npz'], 'Filetype '+filetype+' not understood'
    if filetype in ['.txt', '.dat']:
        idnames = np.array(open(clustermat, 'r').readline().strip().split()[1:])
        idmat = np.genfromtxt(clustermat, dtype = float)
    else:
        sfile = np.load(clustermat)
        idmat = sfile['identmat']
        idnames = sfile['names'].astype(str)
    if splitnames is not None:
        for i in range(len(idnames)):
            idnames[i] = idnames[i].split(splitnames)[namepart]
            
    idmat = 1. - idmat/maxid
    return idnames, idmat


if '--seqidfile' in sys.argv:
    # sequence identity between all RBRs of 53 eukaryotes
    simfile = sys.argv[sys.argv.index('--seqidfile')+1]
    labels, idmat = readinclustermat(simfile, splitnames = '||')
    outname = os.path.splitext(simfile)[0]

if '--rbplist' in sys.argv:
    rlist = np.genfromtxt(sys.argv[sys.argv.index('--rbplist')+1], dtype = str)
    mask = np.isin(labels, rlist)
    labels, idmat = labels[mask], idmat[mask][:,mask]

if '--domainclass' in sys.argv:
    domainfile = np.genfromtxt(sys.argv[sys.argv.index('--domainclass')+1], dtype = str)
    domains = []
    for l, label in enumerate(labels):
        domains.append(domainfile[list(domainfile[:,1]).index(label), -1])
    domains = np.array(domains)
    
    if '--domain_sensitive' in sys.argv:
        for d, domain in enumerate(domains):
            for e in range(d+1,len(domains)):
                cd = True
                if 'RRM' in domain and 'RRM' in domains[e]:
                    cd = False
                elif 'KH' in domain and 'KH' in domains[e]:
                    cd = False
                # set distances to 1. if two proteins don't share any domain types
                if cd:
                    idmat[d,e] = idmat[e,d] = 1.
    
    elif '--domain_select' in sys.argv:
        # remove all sequences that contain a different domain type
        selecteddomain = sys.argv[sys.argv.index('--domain_select')+1]
        outname += 'selectrbd'+selecteddomain
        protmask = np.ones(len(idmat)) == 1
        for d, domain in enumerate(domains):
            if selecteddomain not in domain:
                protmask[d] = False
        labels, idmat, domains = labels[protmask], idmat[:, protmask][protmask], domains[protmask]
    
     
    if '--domainadjust' in sys.argv:
        domaincounts = np.array([len(dc) for dc in domains ])
        domaincountmat = np.ones((len(domaincounts), len(domaincounts))) * domaincounts
        domaindiffmat = np.absolute(domaincountmat - domaincountmat.T) <= 1
        domaincountmat = np.amin(np.array([domaincountmat, domaincountmat.T]),axis = 0)/np.amax(np.array([domaincountmat, domaincountmat.T]),axis = 0)
        domaincountmat[domaindiffmat] = 1.
        idmat = 1.-(idmat-1.)*domaincountmat
    
    
if '--cutoff' in sys.argv:
    cutoff = float(sys.argv[sys.argv.index('--cutoff')+1])
else: 
    cutoff = 0.

def groups_inlinkage(linkmatrix, cut):
    linksets = [[i] for i in range(len(linkmatrix) +1)]
    for l, link in enumerate(linkage_matrix):
        linksets.append(np.append(linksets[int(link[0])], linksets[int(link[1])]))
    sets = linksets[len(linkmatrix) + 1:]
    groups = []
    isin = []
    for l, link in enumerate(linkage_matrix):
        if link[2] > cut and not np.isin(sets[l], isin).any():
            groups.append(sets[l])
            isin = np.append(isin, sets[l])
    return groups
    
# choose the agglomerative clustering algorithm that builds the tree
if '--linkage' in sys.argv:
    linkagemeth = sys.argv[sys.argv.index('--linkage')+1]
else:
    linkagemeth = 'single'


if linkagemeth == 'nj':
    print('build nj tree')
    dm = DistanceMatrix(idmat, labels)
    newick_str = nj(dm, result_constructor=str)
    tree = Phylo.read(StringIO(newick_str), "newick")
    testgroups = []
    isin = []
    for clade in tree.depths():
        if tree.depths()[clade] > cutoff:
            ngroup = [leaf.name for leaf in clade.get_terminals()]
            if not np.isin(ngroup, isin).any() and len(ngroup) > 1:
                testgroups.append(np.where(np.isin(labels, ngroup))[0])
                isin = np.append(isin,ngroup)
    
    
# Use precalculated tree (for example from ClustalW)
elif linkagemeth == 'precalc':
    tree = Phylo.read(sys.argv[sys.argv.index('--linkage')+2], "newick")
    leafs = tree.get_terminals()
    sorting = []
    for leaf in leafs:
        if '__' in leaf.name:
            leaf.name = leaf.name.split('__')[0]
        sorting.append(leaf.name)
    sorting = np.isin(labels, sorting)
    
    classmat = classmat[sorting]
    classsize = classsize[sorting]
    confidence = confidence[sorting]
    evolution = evolution[sorting]
    parent = parent[sorting]
    
    jpledistance = jpledistance[sorting]
    latentnames = latentnames[sorting]
    labels = labels[sorting]
    idmat = idmat[sorting][:,sorting]
    labels = labels[sorting]
    domains = domains[sorting]
    bestid = bestid[sorting]
        
else:    
    linkage_matrix = linkage(idmat[np.triu_indices(len(idmat),1)], method=linkagemeth)
    testgroups = groups_inlinkage(linkage_matrix, cutoff)

obj = open(outname+'_'+linkagemeth+str(cutoff)+'_sets.txt', 'w')
print('Writing', outname+'_'+linkagemeth+str(cutoff)+'_sets.txt')
for t, test in enumerate(testgroups):
    print( 'Testset' , t, len(test))
    train = np.delete(np.arange(len(labels), dtype = int), test)
    obj.write('###Set '+str(t)+'\n'+'##Train:'+'\n'+' '.join(labels[train])+'\n'+' '.join(labels[train])+'\n'+'##Test:'+'\n'+' '.join(labels[test])+'\n'+' '.join(labels[test])+'\n')
maxtotrain = []
for t, test in enumerate(testgroups):
    train = np.delete(np.arange(len(labels), dtype = int), test)
    for te in test:
        print(labels[te], int(np.amax(1.-idmat[te][train]) * 100.))
        maxtotrain.append(int(np.amax(1.-idmat[te][train]) * 100.))

figmaxid = plt.figure(figsize = (4,4))
axmax = figmaxid.add_subplot(111)
axmax2 = axmax.twinx()
axmax.hist(maxtotrain, bins = 20, label = 'Histogram', alpha = 0.6)
axmax2.hist(maxtotrain, bins = 20, density=True, histtype='step', cumulative=True, color = 'k', lw = 2., label = 'Cummulative')
axmax.spines['top'].set_visible(False)
axmax.spines['right'].set_visible(False)
axmax2.spines['top'].set_visible(False)
axmax2.spines['left'].set_visible(False)
#axmax2.spines['bottom'].set_visible(False)
axmax.set_xlabel('Max SeqID to train')
axmax2.set_ylabel("Percent test RBPs")
axmax2.set_yticks([0,0.25, 0.5, 0.75, 1.])
axmax2.set_yticklabels([0,25, 50, 75, 100])
axmax2.grid(ls = '--')
axmax.grid(axis = 'x', ls = '--')
axmax.set_ylabel('Number test RBPs')
   
fig = plt.figure(figsize = (6,0.2*len(labels)), dpi = 20)
axden = fig.add_subplot(111)
axden.spines['top'].set_visible(False)
axden.spines['left'].set_visible(False)
axden.spines['right'].set_visible(False)
axden.tick_params(which = 'both', left = False, labelleft = False, labelright = True)


if '--protnamefile' in sys.argv:
    pnames = np.genfromtxt(sys.argv[sys.argv.index('--protnamefile')+1], dtype = str, delimiter = ' ', )
    protnames = []
    for l, lab in enumerate(labels):
        protnames.append(pnames[list(pnames[:,1]).index(lab), 2]+"__"+domains[l])
    protnames = np.array(protnames)
else:
    protnames = labels
    
    
with plt.rc_context({'lines.linewidth': 1.5}):
    if linkagemeth == 'nj' or linkagemeth == 'precalc':
        leafs = tree.get_terminals()
        lengths = []
        for leaf in leafs:
            lengths.append(tree.distance(leaf))
        print (np.amax(lengths))
        
        def lfunc(cl):
            return None
        Phylo.draw(tree, label_func=lfunc, axes = axden, do_show = False)
        sortree = tree.get_terminals('postorder')
        presortree = labels.tolist() #np.arange(len(sortree),dtype = int).astype(str)
        sorting = []
        for so in sortree:
            if '__' in so.name:
                so.name = so.name.split('__')[0]
            sorting.append(presortree.index(so.name))
        axden.set_ylim([0.5,len(labels)+0.5])
        axden.set_xlim([axden.get_xlim()[0], np.amax(lengths)])
        axden.set_yticks(np.sort(sorting)+1)
        axden.set_yticklabels(protnames[sorting])
        axden.plot([cutoff,cutoff],[0,len(labels)*10], color = 'r', lw = 1.)
    else:
        dn = dendrogram(linkage_matrix, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = axden)
        sorting = dn['leaves']
        axden.set_ylim([0,len(labels)*10])
        axden.set_xticks([0.3,0.5,0.7,1.])
        axden.set_xticklabels(['70%','50%','30%','0%'])
        axden.set_xlabel('Sequence ID')
        axden.grid(axis ='x')
        axden.plot([cutoff,cutoff],[0,len(labels)*10], color = 'r', lw = 1.)
        axden.set_yticks((np.sort(sorting)+0.5)*10)
        axden.set_yticklabels(protnames[sorting])

dpi = 100
if '--dpi' in sys.argv:
    dpi = int(sys.argv[sys.argv.index('--dpi')+1])
if '--savefig' in sys.argv:
    print( outname+'_tree_'+linkagemeth+str(cutoff)+'.svg')
    fig.savefig(outname+'_tree_'+linkagemeth+str(cutoff)+'.svg',dpi = dpi, bbox_inches = 'tight')
    figmaxid.savefig(outname+'_identitydist'+linkagemeth+str(cutoff)+'.svg',dpi = 200, bbox_inches = 'tight')
else:
    plt.show()



