import numpy as np
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, wilcoxon
from sklearn.linear_model import LinearRegression
from scipy.stats import ttest_ind
from scipy.spatial.distance import cdist, cosine
from scipy import stats
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib import cm
from scipy.cluster.hierarchy import dendrogram, linkage

   
## remove branches and leafs that are not in prots
def cleanz(Z, cutoff, prots):
    orleaf = np.arange(len(prots))
    branches = list(np.arange(len(prots), len(prots)+len(Z)))
    incluster = Z[:,2]  < cutoff
    clusters = []
    clusterrep = []
    combbranch = []
    while True:
        start = np.where(incluster)[0][-1]
        cluster = []
        i = start
        ci = branches[i]
        clusterrep.append(ci)
        branchloc = [i]
        branchcover = [0]
        while True:
            if Z[i,0] in cluster:
                cluster.append(Z[i,1])
                branchcover[branchloc.index(i)] += 1
                if Z[i,1] >= len(prots):
                    ci = int(Z[i,1])
                    i = branches.index(ci)
                    branchcover.append(0)
                    branchloc.append(i)
                else:
                    bp = np.where(np.array(branchcover) == 0)[0]
                    if len(bp) == 0:
                        break
                    else:
                        i = branchloc[bp[-1]]
            else:
                cluster.append(Z[i,0])
                if Z[i,0] >= len(prots):
                    ci = int(Z[i,0])
                    i = branches.index(ci)
                    branchcover.append(0)
                    branchloc.append(i)
        cluster = np.array(cluster)
        clusters.append(cluster[cluster < len(prots)].astype(int))
        incluster[branchloc] = False
        combbranch.append(branchloc)
        if np.sum(incluster) == 0:
            break

    assigned = np.concatenate(clusters)
    for o in orleaf:
        if o not in assigned:
            clusters.append([o])
            clusterrep.append(o)

    keepbranch = np.delete(np.arange(len(Z)), np.concatenate(combbranch))
    branchtrans = []
    bc = len(clusters)
    for j in range(len(Z)):
        if j in keepbranch:
            branchtrans.append([bc, branches[j]])
            bc += 1
    znew = Z[keepbranch]
    zred = np.copy(znew)
    for r, rep in enumerate(clusterrep):
        znew[:, :2][zred[:, :2] == rep] = r

    for bt in branchtrans:
        znew[:, :2][zred[:, :2] == bt[1]] = bt[0]

    return clusters, znew



if __name__ == '__main__':
    expressionfile = np.load(sys.argv[1])
    
    expression = expressionfile['pkm']
    tissues = expressionfile['tissues']
    genes = expressionfile['genes']

    valuename = sys.argv[3]
    
    tfile = open(sys.argv[2], 'r').readlines()
    tissueid = []
    for l, line in enumerate(tfile):
        if l> 0:
            line = line.split('\t')
            tissueid.append([line[8], line[10]])
    tissueid = np.array(tissueid)
    uniquetissue, tn = np.unique(tissueid[:, 1], return_index = True)
    tissueloc = []
    
    for u, utiss in enumerate(uniquetissue):
        tissueids = tissueid[tissueid[:,1] == utiss, 0]
        loc = np.where(np.isin(tissues, tissueids))[0]
        tissueloc.append(loc)
    tissueloc = np.array(tissueloc)
    
    expression1 = expression[:, tissueloc.T[0]]
    expression2 = expression[:, tissueloc.T[1]]
    
    expression = (expression1 + expression2)/2.
    
    rbpcorrelation = []
    for g in range(len(genes)):
        rbpcorrelation.append(1.-cosine(expression1[g], expression2[g]))
      
    rbpcorrelation = np.array(rbpcorrelation)
    
    tissuecorrelation = 1. - cdist(expression1.T, expression2.T, 'correlation')
    selftissuecorrelation = 1. - cdist(expression.T, expression.T, 'correlation')
    
 
    
    # plot tissue similarity with matrix and dendogram 
    bottom = 0.1
    height = 0.65
    top = bottom+height
    heightphylo = 0.15
    left = 0.1
    width = 0.8
    
    vmin = 0.
    vmax = 1.
    
    cmapgs = np.concatenate([[[1.,1.,1.,1.] for i in range(20)] , [cm.Oranges(0.1) for i in range(20)], [cm.Oranges(0.3) for i in range(20)], [cm.Oranges(0.6) for i in range(20)], [cm.Oranges(1.) for i in range(20)]], axis = 0)
    colmap = ListedColormap(cmapgs)
        
    fig = plt.figure(figsize = (12,16), dpi = 50)
    ac = fig.add_subplot(222)
    ac.set_position([left+width+0.025, top+0.5*heightphylo, 0.025, 0.5*heightphylo])
    ac.imshow(np.linspace(vmin,vmax,100).reshape(-1,1), origin = 'lower', cmap = colmap, aspect = 'auto', vmin = vmin, vmax = vmax)
    ac.tick_params(left = False, labelleft = False, right = True, labelright = True, labelbottom = False, bottom = False, top = False, labeltop = False)
    ac.set_yticks([0,50,99])
    ac.set_yticklabels([vmin, (vmin + vmax)/2., vmax])
    ac.spines['top'].set_visible(False)
    ac.spines['left'].set_visible(False)
    ac.spines['right'].set_visible(False)
    ac.spines['bottom'].set_visible(False)        
    ac.set_title('Pearson R')
    
    ax0 = fig.add_subplot(121)
    ax0.set_position([left,top,width,heightphylo])
    zcore = cdist(expression.T, expression.T, 'correlation')
    Z = linkage(zcore[np.triu_indices(len(zcore),1)], 'single')
    
    corrcut = 0.87
    mclustercomp, Zred = cleanz(Z, 1.-corrcut, uniquetissue)
    for m, mcomp in enumerate(mclustercomp):
        for mc in mcomp:
            print uniquetissue[mc]+'\t'+str(m)
        
            
    dn = dendrogram(Z, orientation = 'top', color_threshold=1.-corrcut, above_threshold_color='k', ax = ax0)
    limx = ax0.get_xlim()
    ax0.plot([np.amin(limx),np.amax(limx)], [1.-corrcut,1.-corrcut], c = 'r', ls = '--')
    sorting = dn['leaves']
    ax0.tick_params(left = True, labelleft = True, right = False, labelright = False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.set_ylabel("Pearson")
    ax0.set_yticks([0,1.-corrcut,0.4])
    ax0.set_yticklabels([1,corrcut,0.6])
    
    ax = fig.add_subplot(122)
    ax.set_position([left,bottom,width,height])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(which = 'both', left = True, labelbottom = False, bottom = False, labelleft = True, labelright = False)
    im = ax.imshow(selftissuecorrelation[sorting][:, sorting], origin = 'lower', cmap = colmap, aspect = 'auto', vmin = vmin, vmax = vmax)
    ax.set_yticks(np.arange(len(uniquetissue)+1)-0.5, minor = True)
    ax.set_xticks( np.arange(len(uniquetissue)+1)-0.5, minor = True)
    ax.grid(color='silver', linewidth=0.5, which = 'minor')
    #ax.set_xticks(np.arange(len(uniquetissue)))
    #ax.set_xticklabels(uniquetissue[sorting[::-1]], rotation = 60, ha = 'right', va = 'top' ) 
    ax.set_yticks(np.arange(len(uniquetissue)))
    ax.set_yticklabels(uniquetissue[sorting]) 
    
    sorting = sorting
    
    # plot distribution of expression for each tissue
    figtdist = plt.figure(figsize = (4., len(uniquetissue)*0.25), dpi = 50)
    axtd = figtdist.add_subplot(111)
    axtd.spines['top'].set_visible(False)
    #axtd.spines['bottom'].set_visible(False)
    axtd.spines['right'].set_visible(False)
    axtd.set_xlabel(valuename)
    axtd.boxplot(list(expression.T[sorting]), showfliers = False, vert = False)
    ylimd = axtd.get_ylim()
    axtd.plot([np.mean(expression),np.mean(expression)],ylimd, c = 'grey')
    
    axtd.set_yticks(np.arange(1, len(uniquetissue)+1))
    axtd.set_yticklabels(uniquetissue[sorting])
    #axtd.set_xlim([0,0.025])
    
    
    # plot boxplot for row and colum each 
    # plot where self correlation lies
    # sort by selfcorrelation # 
    figtdists = plt.figure(figsize = (4, len(uniquetissue)*0.25), dpi = 50)
    axtds = figtdists.add_subplot(111)
    axtds.spines['top'].set_visible(False)
    #axtd.spines['bottom'].set_visible(False)
    axtds.spines['right'].set_visible(False)
    tiscor =tissuecorrelation[sorting][:, sorting]
    np.fill_diagonal(tiscor, 0.5)
    b1 = axtds.boxplot(list(tiscor), positions = np.arange(len(uniquetissue)), widths = 0.35, patch_artist = True, vert = False)
    b2 = axtds.boxplot(list(tiscor.T), positions = np.arange(len(uniquetissue))+0.35, widths = 0.35, patch_artist = True, vert = False)
    for patch in b1['boxes']:
        patch.set(facecolor='tan')
    for patch in b2['boxes']:
        patch.set(facecolor='powderblue')
    axtds.scatter(tissuecorrelation.diagonal()[sorting], np.arange(len(uniquetissue))+0.175, c = 'r') #, label = 'Tissue replicate')
    axtds.set_yticks(np.arange(0.175, len(uniquetissue)+0.175))
    axtds.set_yticklabels(uniquetissue[sorting])
    axtds.set_xlabel('Pearson replicates')
    #axtds.legend()
    
    
    # plot distribution of correlation for genes
        # color by median expression
    figpdist = plt.figure(figsize = (4, 2.5), dpi = 300)
    axpd = figpdist.add_subplot(111)
    axpd.spines['top'].set_visible(False)
    #axpd.spines['left'].set_visible(False)
    axpd.spines['right'].set_visible(False)
    axpd.tick_params(which = 'both')
    bins = np.linspace(np.amin(rbpcorrelation), np.amax(rbpcorrelation), 33)
    axpd.spines['left'].set_color('goldenrod')
    medianexpression = np.median(expression, axis = 1)
    bincenter = []
    medianexp = []
    for b in range(len(bins)-1):
        bincenter.append((bins[b]+bins[b+1])/2.)
        medianexp.append(np.mean(medianexpression[(rbpcorrelation > bins[b])*(rbpcorrelation < bins[b+1])]))
    medianexp = np.nan_to_num(medianexp)
    signy = 1.
    if np.array(medianexp)[np.argmax(np.absolute(medianexp))] <0:
        signy = -1.
    sns.distplot(rbpcorrelation, bins = bins, ax = axpd, label = 'Replicate correlation (genes)')
    limy = axpd.get_ylim()
    
    axpd.set_yticks([0,limy[1]*0.9/2., limy[1]*0.9])
    medianexp = np.around(medianexp, 4)
    medianmask = ~np.isnan(medianexp)
    medianexp = medianexp[medianmask]
    bincenter = np.array(bincenter)[medianmask]
    axpd.set_yticklabels(np.around([signy*np.amin(signy*medianexp), signy*(np.amin(signy*medianexp)+np.amax(signy*medianexp))/2.,signy* np.amax(signy*medianexp)],4))
    medianexp = limy[1]*0.9*(medianexp-signy*np.amin(signy*medianexp))/(signy*np.amax(signy*medianexp)-signy*np.amin(signy*medianexp))
    axpd.scatter(bincenter, medianexp, c = 'goldenrod', label = 'Mean rpkm per bin')
    axpd.set_xlabel('Pearson genes (replicates)')
    axpd.set_ylabel('Mean '+valuename)
    axpd.legend(prop={'size':8})
    
    

    

    # plot difference distribution and sort wilcoxon, maybe color by wilcoxon
    diffexpression = expression1 - expression2
    #diffexpression = diffexpression/medianexpression[:, None]
    meandiff = np.mean(diffexpression, axis = 0)
    diffexpression = np.sign(meandiff)[None, :]*diffexpression
    figtdif = plt.figure(figsize = (4, len(uniquetissue)*0.25))
    axtdf = figtdif.add_subplot(111)
    axtdf.spines['top'].set_visible(False)
    #axtd.spines['bottom'].set_visible(False)
    axtdf.spines['right'].set_visible(False)
    axtdf.set_xlabel('Delta genes '+valuename)
    axtdf.boxplot(list(diffexpression.T[sorting]), vert = False)
    axtdf.set_yticks(np.arange(1, len(uniquetissue)+1))
    axtdf.plot( [0,0],[0,len(uniquetissue)+1], color = 'grey')
    axtdf.set_yticklabels(uniquetissue[sorting])
    
    
    

    
    
    
    
    fig.savefig('Arabidopsis_tissue_expression_correlation_'+valuename+'.jpg', dpi = 250, bbox_inches = 'tight')
    figtdist.savefig('Arabidopsis_tissue_expression_distributionzoom_'+valuename+'.jpg', dpi = 250, bbox_inches = 'tight')
    figpdist.savefig('Arabidopsis_gene_expression_correlation_distribution_'+valuename+'.jpg', dpi = 250, bbox_inches = 'tight')
    figtdists.savefig('Arabidopsis_tissue_replicate_correlation_'+valuename+'.jpg', dpi = 250, bbox_inches = 'tight')
    figtdif.savefig('Arabidopsis_tissue_replicatedifference_correlation_'+valuename+'.jpg', dpi = 250, bbox_inches = 'tight')
    
    
    
    
    #plt.show()
    
    
    ### stats for 3'utrs
        # len vs n binding sites
        # len vs n rbp




    
    
