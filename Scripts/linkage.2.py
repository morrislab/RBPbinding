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

# make heatmap with sequence id (overlap of p-vectors) and jple distance to other species, include all species on tree.
# make barplot for each species with different coloring for number of species the RBP clusters are in. See which bar is extraordonary big


# esiest solution, use different seq id measure for tree
# biopython has neighbor joining method
# To do:  use biopython to readin clustalomega identity matrix or to run in python
# Use sequence identity of similarity and create own NJ algorithm to create linkage matrix
# somehow plot with ete3 or biopython (don't know who to do yet)
# helpful tools:
#https://www.biotite-python.org/apidoc/biotite.sequence.phylo.neighbor_joining.html
#https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html
#https://github.com/sness/courses/blob/master/neighbour-joining.py
#http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html
#https://en.wikipedia.org/wiki/Neighbor_joining
#https://biopython.org/wiki/Phylo

# node1, node2, distance, number of nodes    
#print(len(idmat), linkage_matrix[:30])
'''
#for neighbor joining
from skbio import DistanceMatrix
from skbio.tree import nj
from Bio import Phylo
from io import StringIO

dm = DistanceMatrix(data)
newick_str = nj(dm, result_constructor=str)
tree = Phylo.read(StringIO(newick_str), "newick")
fig = plt.figure()
ax = fig.add_subplot(111)
def lfunc(cl):
    return None
Phylo.draw(tree, label_func=lfunc, axes = ax)
# get order of tree
sort = tree.get_terminals('postorder')
sorting = []
for so in sort:
    sorting.append(labels.index(so.name))

'''


# Read in different files
def readinclustermat(clustermat, twonames = False, maxid = 100.):
    sfile = np.load(clustermat)
    idmat = sfile['identmat']
    idnames = sfile['names'].astype(str)
    idmat = 1. - idmat/maxid
    if twonames:
        idnames2 = sfile['names2'].astype(str)
        return idnames, idmat, idnames2
    return idnames, idmat






if '--seqidfile' in sys.argv:
    # sequence identity between all RBRs of 53 eukaryotes
    simfile = sys.argv[sys.argv.index('--seqidfile')+1]
    labels, idmat = readinclustermat(simfile)
    if '--outname' in sys.argv:
        outname = sys.argv[sys.argv.index('--outname')+1]
    else:
        outname = os.path.splitext(simfile)[0]
    
    protnames = []
    domains = []
    
    for l, label in enumerate(labels):
        protnames.append(label.split('__')[0])
        domains.append(label.split('__')[3].split('-'))
        # domainplot adds a column for the representation of the domain composition
        if '--domainplot' in sys.argv:
            labels[l] = label.split('__')[2]
        else:
            labels[l] = label.split('__')[2]+'_'+label.split('__')[3]
    protnames = np.array(protnames)
    domains = np.array(domains)

    # Names of species that are compared
    simnames = np.array(sys.argv[sys.argv.index('--seqidfile')+3].split(',')) # arabi_t,arabi_l,briga    
    # latentfiles contain the latent representations for the RBRs in all compared species
    latentfiles = sys.argv[sys.argv.index('--seqidfile')+2].split(',') # need to be in same order as simfiles
    
    latent = []
    latentnames = []
    for l, latfile in enumerate(latentfiles):
        slat = np.load(latfile)
        slatnam = slat['names'].astype(str)
        slatsort = np.argsort(slatnam)
        latent.append(slat['profiles'][slatsort])
        latentnames.append(slatnam[slatsort])
    
    # sort protnames and only keep the ones that are also in latentnames
    protmask = np.argsort(protnames)
    protmask = protmask[np.isin(protnames[protmask], np.concatenate(latentnames))]
    labels, idmat, domains, protnames = labels[protmask], idmat[:, protmask][protmask], domains[protmask], protnames[protmask]
    
    if '--domain_sensitive' in sys.argv:
        dclass = np.zeros((2,len(domains)))
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
        labels, idmat, domains, protnames = labels[protmask], idmat[:, protmask][protmask], domains[protmask], protnames[protmask]
    
    
    jpledistance = []
    # iterate over species that were included
    for l, lat in enumerate(latent):
        closedist = []
        for k, lat2 in enumerate(latent):
            # compute distances between all RBPs of two species
            cist = cdist(lat, lat2, 'cosine')
            # only include the rbps that are also in protnames
            cist = cist[np.isin(latentnames[l], protnames)][:, np.isin(latentnames[k], protnames)]
            # determine the sequence identities between the proteins of the two species
            cidmat = idmat[np.isin(protnames,latentnames[l])][:, np.isin(protnames,latentnames[k])]
            # RRM and KH domains all have 1. in idmat if domain_sensitive, # we set the cosine distance to maximum between KH and RRM domains as well
            cist[cidmat>0.9] = 2.
            # within the same species, we measure the distance to the closest RBP in the species
            if l == k:
                closedist.append(np.sort(cist, axis = 1)[:,1])                
            else:
                closedist.append(np.amin(cist, axis = 1))
        jpledistance.append(closedist)
    jpledistance = np.concatenate(jpledistance, axis = 1).T
    
    # Find the closest protein sequence in the other species
    diffspecies = np.zeros(len(labels), dtype = int)
    bestid = np.zeros((len(labels), len(latentnames)))
    for i in range(len(latentnames)):
        maski = np.isin(protnames, latentnames[i])
        diffspecies[maski] = i
        for j in range(len(latentnames)):
            maskj = np.isin(protnames, latentnames[j])
            if i != j:
                bestid[maski,j] = np.amin(idmat[maski][:,maskj], axis = 1)
            else:
                bestid[maski,j] = np.sort(idmat[maski][:,maskj], axis = 1)[:,1]
    bestid = (1.-bestid)*100.
    
    # Create short names the species that are tested
    shortspecies = []
    for sn, sname in enumerate(simnames):
        if '.' in sname:
            shortspecies.append(sname.split('.')[0][0]+ sname.split('.')[1][0])
        elif '_' in sname:
            shortspecies.append(sname.split('_')[0][0]+ sname.split('_')[1][0])
    outname += '-'.join(shortspecies)
    
    # sort latentnames after protnames
    latentnames = np.concatenate(latentnames)
    latentnames = latentnames[np.isin(latentnames, protnames)]
    sort = []
    for label in protnames:
        sort.append(list(latentnames).index(label))
    jpledistance = jpledistance[sort]
    
    





    
if '--cluster' in sys.argv:
    clusterfile = np.genfromtxt(sys.argv[sys.argv.index('--cluster')+1], dtype = str)
    species = np.array(open(sys.argv[sys.argv.index('--cluster')+1], 'r').readline().strip().split()[14:])
    clusters = []
    sort = np.ones(len(protnames)) == 1
    for p, prot in enumerate(protnames):
        if prot in clusterfile[:,0]:
            clusters.append(clusterfile[list(clusterfile[:,0]).index(prot),2:])
        else:
            sort[p] = False
        #print( maxid[p], clusters[-1])
    clusters = np.array(clusters)
    classmat = clusters[:,[0,1]].astype(int)
    classage = clusters[:,5].astype(float)
    if '--classage' in sys.argv:
        classsize = classage
    else:
        classsize = clusters[:,4].astype(float)
    confidence = np.around(clusters[:,2].astype(float),2)
    evolution = clusters[:,9].astype(float)
    
    jpledistance = jpledistance[sort]
    latentnames = latentnames[sort]
    protnames = protnames[sort]
    idmat = idmat[sort][:,sort]
    labels = labels[sort]
    domains = domains[sort]
    bestid = bestid[sort]
    

 
# choose the agglomerative clustering algorithm that builds the tree
if '--linkage' in sys.argv:
    linkagemeth = sys.argv[sys.argv.index('--linkage')+1]
else:
    linkagemeth = 'single'


if linkagemeth == 'nj':
    print('build nj tree')
    dm = DistanceMatrix(idmat, protnames)
    newick_str = nj(dm, result_constructor=str)
    tree = Phylo.read(StringIO(newick_str), "newick")
    print('done')
# Use precalculated tree (for example from ClustalW)
elif linkagemeth == 'precalc':
    tree = Phylo.read(sys.argv[sys.argv.index('--linkage')+2], "newick")
else:    
    linkage_matrix = linkage(idmat[np.triu_indices(len(idmat),1)], method=linkagemeth)

 
 

domaintypes = ['KH', 'RRM']
domaincolor = ListedColormap(['white','firebrick', 'purple'])

mylen = np.vectorize(len)

if '--domainplot' in sys.argv:
    domainmat = np.zeros((len(domains), np.amax(mylen(domains))))
    for d, domain in enumerate(domains):
        for e, dtt in enumerate(domain):
            if dtt in domaintypes:
                domainmat[d,e] = domaintypes.index(dtt)+1

cmapgs = np.concatenate([[[1.,1.,1.] for i in range(3)], sns.color_palette('YlOrBr', 4), [[0.39529412, 0.09582468, 0.] for i in range(3)]], axis = 0)
cmapgs = ListedColormap(cmapgs)
confmap = ListedColormap(np.concatenate([sns.color_palette('Blues_r', 6)[1:],[[1.,1.,1.]]], axis = 0) )
    
fig = plt.figure(figsize = (8,0.25*len(labels)), dpi = 20)
axden = fig.add_subplot(191)
axden.spines['top'].set_visible(False)
axden.spines['left'].set_visible(False)
axden.spines['right'].set_visible(False)
axden.tick_params(which = 'both', left = False, labelleft = False)

with plt.rc_context({'lines.linewidth': 3.}):
    if linkagemeth == 'nj':
        leafs = tree.get_terminals()
        lengths = []
        for leaf in leafs:
            lengths.append(tree.distance(leaf))
        print (np.amax(lengths))
        
        def lfunc(cl):
            return None
        Phylo.draw(tree, label_func=lfunc, axes = axden, do_show = False)
        sortree = tree.get_terminals('postorder')
        presortree = protnames.tolist() #np.arange(len(sortree),dtype = int).astype(str)
        sorting = []
        for so in sortree:
            sorting.append(presortree.index(so.name))
        axden.set_ylim([0.5,len(labels)+0.5])
        print(axden.get_xlim())
        axden.set_xlim([axden.get_xlim()[0], np.amax(lengths)])
    else:
        dn = dendrogram(linkage_matrix, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = axden)
        sorting = dn['leaves']
        axden.set_ylim([0,len(labels)*10])
        axden.set_xticks([0.3,0.5,0.7,1.])
        axden.set_xticklabels(['70%','50%','30%','0%'])
        axden.set_xlabel('Sequence ID')
        axden.grid(axis ='x')
        axden.plot([0.7,0.7],[0,len(labels)*10], color = 'r', lw = 2.)
        axden.plot([0.6,0.6],[0,len(labels)*10], color = 'r', lw = 1., ls = '--')

if '--cluster' in sys.argv:
    axcl = fig.add_subplot(192)
    axm = fig.add_subplot(193)
    axb = fig.add_subplot(195)
    axe = fig.add_subplot(198)
    axco = fig.add_subplot(197)
    
    axden.set_position([0.1,0.1,0.4,0.8])
    axcl.set_position([0.65,0.1,0.05,0.8])
    axe.set_position([0.55,0.1,0.05,0.8])
    axco.set_position([0.6,0.1,0.05,0.8])
    axm.set_position([0.5,0.1,0.05,0.8])
    
    of = 0
    if '--speciesmat' in sys.argv:
        axb.set_position([0.7,0.1,0.6,0.8])
        of = 0.6
    elif '--numbermatrix' in sys.argv:
        axb.set_position([0.7,0.1,0.05,0.8])
        of = 0.05
    else:
        axb.set_position([0.7,0.1,0.3,0.8])
        of = 0.3
    
    of2 = 0.
    of3 = 0.
    if '--seqidfile' in sys.argv:
        axc = fig.add_subplot(196)
        of2 = 0.04*len(shortspecies)
        axc.set_position([0.7+of,0.1,of2,0.8])
        of3 = 0.05*len(shortspecies)
        axc2 = fig.add_subplot(199)
        axc2.set_position([0.7+of+of2,0.1,of3,0.8])
        
        axcb = fig.add_subplot(296)
        axcb.set_position([0.7+of,0.1-0.8/float(len(sorting)),of2,0.8/float(len(sorting))])
        axc2b = fig.add_subplot(299)
        axc2b.set_position([0.7+of+of2,0.1-0.8/float(len(sorting)),of3,0.8/float(len(sorting))])
        
    
    if '--domainplot' in sys.argv:
        axd = fig.add_subplot(194)
        axd.set_position([0.7+of+of2+of3+0.01,0.1,0.18,0.8])
        
    
    
    

    classmat = classmat[sorting]
    oc = -1
    ci = 0
    classchange = []
    for c, cl in enumerate(classmat):
        if cl[0] != oc:
            oc = np.copy(cl[0])
            ci += 1
        classchange.append(ci%2)
    classchange = np.array(classchange).reshape(-1,1)
    
    axm.imshow(classchange,origin = 'lower', cmap = 'Greys', aspect = 'auto', vmin = 0, vmax = 3)
    axm.spines['top'].set_visible(False)
    axm.spines['bottom'].set_visible(False)
    axm.spines['left'].set_visible(False)
    axm.spines['right'].set_visible(False)
    axm.tick_params(which = 'both', left = False, bottom = True, labelleft = False, labelbottom = True)
    axm.set_ylim([-0.5,len(labels)-0.5])
    axm.set_xticks([0])
    axm.set_xticklabels(['Cluster'], rotation = 30, ha = 'right', va = 'top')
    for i in range(len(clusters)):
        axm.text(0,i,str(int(classmat[i][0])), ha = 'center', va = 'center')
    
    evolution = evolution[sorting].reshape(-1,1) 
    axe.imshow(evolution, origin = 'lower', cmap = ListedColormap(cm.tab20([13,12,10,6,8])), aspect = 'auto')
    axe.spines['top'].set_visible(False)
    axe.spines['bottom'].set_visible(False)
    axe.spines['left'].set_visible(False)
    axe.spines['right'].set_visible(False)
    axe.tick_params(which = 'both', left = False, bottom = True, labelleft = False, labelbottom = True)
    axe.set_ylim([-0.5,len(labels)-0.5])
    axe.set_xticks([0])
    axe.set_xticklabels(['Evolution'], rotation = 30, ha = 'right', va = 'top')
    for i in range(len(evolution)):
        axe.text(0,i,str(int(evolution[i][0])), ha = 'center', va = 'center')
    
    confidence = confidence[sorting].reshape(-1, 1)
    axco.imshow(confidence ,origin = 'lower', cmap = confmap, aspect = 'auto', vmin = 0., vmax = 0.6)
    axco.spines['top'].set_visible(False)
    axco.spines['bottom'].set_visible(False)
    axco.spines['left'].set_visible(False)
    axco.spines['right'].set_visible(False)
    axco.tick_params(which = 'both', left = False, bottom = True, labelleft = False, labelbottom = True)
    axco.set_ylim([-0.5,len(labels)-0.5])
    axco.set_xticks([0])
    axco.set_xticklabels(['Confidence'], rotation = 30, ha = 'right', va = 'top')
    for i in range(len(clusters)):
        axco.text(0,i,str(confidence[i][0]), ha = 'center', va = 'center')
    
    
    
    if '--seqidfile' in sys.argv:
        diffspecies = diffspecies[sorting].reshape(-1,1)
        axcl.imshow(diffspecies, origin = 'lower', cmap = 'Paired_r', aspect = 'auto', vmin = 0, vmax = 8)
        axcl.set_xticks([0])
        axcl.set_xticklabels(['Species'], rotation = 30, ha = 'right', va = 'top')
        for i, ds in enumerate(diffspecies):
            axcl.text(0,i, str(ds[0]+1), color = 'k', ha = 'center', va = 'center')#shortspecies[ds[0]]
            
    else:
        axcl.imshow(classmat[:,[1]],origin = 'lower', cmap = 'Reds', aspect = 'auto', vmin = 0, vmax = 1)
        axcl.set_xticks([0])
        axcl.set_xticklabels(['Measured'], rotation = 30, ha = 'right', va = 'top')
    axcl.spines['top'].set_visible(False)
    axcl.spines['bottom'].set_visible(False)
    axcl.spines['left'].set_visible(False)
    axcl.spines['right'].set_visible(False)
    axcl.tick_params(which = 'both', left = False, bottom = True, labelleft = False, labelbottom = True)
    axcl.set_ylim([-0.5,len(labels)-0.5])
    
    
    if '--domainplot' in sys.argv:
        axd.imshow(domainmat[sorting],origin = 'lower', cmap = domaincolor, aspect = 'auto', vmin = 0, vmax = 2)
        axd.spines['top'].set_visible(False)
        axd.spines['bottom'].set_visible(False)
        axd.spines['left'].set_visible(False)
        axd.spines['right'].set_visible(False)
        axd.tick_params(which = 'both', left = False, bottom = True, labelleft = False, labelbottom = True)
        axd.set_ylim([-0.5,len(labels)-0.5])
        axd.set_xticks([np.around(len(domainmat[0])/2.,0)])
        axd.set_xticklabels(['Domain class'], rotation = 30, ha = 'right', va = 'top')
        axd.set_xticks(np.arange(len(domainmat[0]))+0.5, minor = True)
        axd.set_yticks(np.arange(len(domainmat))+0.5, minor = True)
        axd.grid(which= 'minor', color = 'white', lw = 3.)
        
    
    if '--speciesmat' in sys.argv:
        axb.imshow(clusters[sorting][:,6:].astype(int),origin = 'lower', cmap = 'Greys', aspect = 'auto', vmin = 0, vmax = 1)
        axb.set_ylim([-0.5,len(labels)-0.5])
        axb.spines['top'].set_visible(False)
        axb.spines['left'].set_visible(False)
        axb.spines['right'].set_visible(False)
        axb.tick_params(which = 'both', left = False, labelleft = False, right = True, labelright = True)
        axb.set_xticks(np.arange(len(species)))
        axb.set_xticklabels(species, rotation = 90)
        
        axb.set_yticks(np.arange(len(labels))+0.5)
        axb.set_yticklabels(labels[sorting])
        
    elif '--numbermatrix':
        classsize = classsize[sorting].reshape(-1,1).astype(int)
        if '--classage' in sys.argv:
            vmin, vmax = 1, 300
        else:
            vmin, vmax = 1, 3
        axb.imshow(classsize,origin = 'lower', cmap = 'Greys', aspect = 'auto', vmin = vmin, vmax = vmax)
        axb.set_ylim([-0.5,len(labels)-0.5])
        axb.spines['top'].set_visible(False)
        axb.spines['left'].set_visible(False)
        axb.spines['right'].set_visible(False)
        axb.tick_params(which = 'both', left = False, labelleft = False, right = True, labelright = True)
        axb.set_xticks([0])
        axb.set_xticklabels(['Number of species\npresent'], rotation = 90)
        for c, clsz in enumerate(classsize):
            if clsz[0] <= 1:
                axb.text(0, c, str(clsz[0]), color = 'grey', va = 'center' , ha = 'center')
            else:
                axb.text(0, c, str(clsz[0]), color = 'white', va = 'center' , ha = 'center')

    else:

        axb.barh(np.arange(len(labels))+0.5, classsize[sorting], color = 'grey', edgecolor = 'k')
        axb.set_ylim([0,len(labels)])
        axb.spines['top'].set_visible(False)
        axb.spines['left'].set_visible(False)
        axb.spines['right'].set_visible(False)
        axb.tick_params(which = 'both', left = False, labelleft = False, right = True, labelright = True)
        #axb.set_xscale('log')
        axb.set_xlabel('Number of species\npresent')
        axb.set_xticks([1,10,25,50])
        axb.set_xticklabels([1,10,25,50])

        
    if '--seqidfile' in sys.argv:
        bestid = bestid[sorting].astype(int)
        axc.imshow(bestid,origin = 'lower', cmap = cmapgs, aspect = 'auto', vmin = 0, vmax = 100)
        axc.set_ylim([-0.5,len(labels)-0.5])
        axc.spines['top'].set_visible(False)
        axc.spines['left'].set_visible(False)
        axc.spines['right'].set_visible(False)
        axc.tick_params(which = 'both', left = False, labelleft = False, right = True, labelright = True)
        
        
        for bi, beid in enumerate(bestid):
            for bj, bid in enumerate(beid):
                axc.text(bj,bi,str(bid), color = 'grey', va = 'center' , ha = 'center')
        
        axcb.imshow(np.arange(len(shortspecies)).reshape(1,-1), origin = 'lower', cmap = 'Paired_r', aspect = 'auto', vmin = 0, vmax = 8)
        axcb.set_xticks(np.arange(len(bestid[0])))
        axcb.set_xticklabels(simnames, rotation = 90)
        axcb.spines['top'].set_visible(False)
        axcb.spines['left'].set_visible(False)
        axcb.spines['right'].set_visible(False)
        axcb.tick_params(which = 'both', left = False, labelleft = False)

        axc2b.imshow(np.arange(len(shortspecies)).reshape(1,-1), origin = 'lower', cmap = 'Paired_r', aspect = 'auto', vmin = 0, vmax = 8)
        axc2b.set_xticks(np.arange(len(bestid[0])))
        axc2b.set_xticklabels(simnames, rotation = 90)
        axc2b.spines['top'].set_visible(False)
        axc2b.spines['left'].set_visible(False)
        axc2b.spines['right'].set_visible(False)
        axc2b.tick_params(which = 'both', left = False, labelleft = False)

        
        jpledistance = np.around(jpledistance[sorting],2)
        axc2.imshow(jpledistance,origin = 'lower', cmap = confmap, aspect = 'auto', vmin = 0, vmax = 0.6)
        axc2.set_ylim([-0.5,len(labels)-0.5])
        axc2.spines['top'].set_visible(False)
        axc2.spines['left'].set_visible(False)
        axc2.spines['right'].set_visible(False)
        axc2.tick_params(which = 'both', left = False, labelleft = False, right = True, labelright = True)
        axc2.set_yticks(np.arange(len(labels)))
        axc2.set_yticklabels(labels[sorting])
        
        for bi, beid in enumerate(jpledistance):
            for bj, bid in enumerate(beid):
                axc2.text(bj,bi,str(bid), color = 'k', va = 'center' , ha = 'center')
        
    else:
        axb.set_yticks(np.arange(len(labels))+0.5)
        axb.set_yticklabels(labels[sorting])
    


else:
    axden.tick_params(which = 'both', right = True, labelright = True)
    axden.set_yticklabels(labels[sorting])

print( outname+'_tree_'+linkagemeth+'.jpg')
dpi = 100
if '--dpi' in sys.argv:
    dpi = int(sys.argv[sys.argv.index('--dpi')+1])
fig.savefig(outname+'_tree_'+linkagemeth+'.jpg',dpi = dpi, bbox_inches = 'tight')
#plt.show()



