import numpy as np
import sys, os
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
# cluster on latent distances
from scipy.optimize import least_squares
from scipy import stats
from statsmodels.stats.multitest import multipletests
# import Z-scores 
# import clusters
# get k-mers for clusters
# random order the clusters and plot and approximate function
# assume that every cluster follwos curve and interpolate to a few new clusters with dotted line and variance. 
# Gaussian processes? 

zscores = np.genfromtxt(sys.argv[1], dtype = float)[:,1:]
znames = open(sys.argv[1],'r').readline().strip().split()[1:]

bindingkmers = []
for z, zname in enumerate(znames):
    #bindingkmers.append(zscores[:,z] >= np.sort(zscores[:, z])[-10])
    #bindingkmers.append(zscores[:,z] >= 4.8)
    pval = 2.*(1 - stats.norm.cdf(zscores[:,z]))
    p_valbind = multipletests(pval, alpha=0.001, method='fdr_bh', is_sorted=False, returnsorted=False)
    if np.sum(p_valbind[0]) > 0:
        zcut = np.amin(zscores[:, z][p_valbind[0]])
    else:
        zcut = np.amax(zscores[:, z])
    #print(z, zcut)
    bindingkmers.append(zscores[:,z] >= zcut)
bindingkmers = np.array(bindingkmers)

latent = np.load(sys.argv[2])
lnames = latent['names'].astype(str)
latentrep = latent['profiles']

sort = []
for l, lname in enumerate(lnames):
    sort.append(znames.index(lname))
bindingkmers = bindingkmers[sort]

add = ''
if '--domain_only' in sys.argv:
    domainfile = np.genfromtxt(sys.argv[sys.argv.index('--domain_only')+1], dtype = str, delimiter = ' ')
    domaintype = sys.argv[sys.argv.index('--domain_only')+2]
    add += '-domain'+domaintype
    keep = []
    for l, lname in enumerate(lnames):
        domain = domainfile[list(domainfile[:,1]).index(lname),-1]
        if domaintype in domain:
            keep.append(l)
    latentrep = latentrep[keep]
    lnames = lnames[keep]
    bindingkmers = bindingkmers[keep]
    domaintype = [domaintype]

if '--domain_split' in sys.argv:
    domainfile = np.genfromtxt(sys.argv[sys.argv.index('--domain_split')+1], dtype = str, delimiter = ' ')
    domaintype = sys.argv[sys.argv.index('--domain_split')+2].split(',')
    add += '-domain'+'-'.join(np.array(domaintype))
    keep = []
    dtype = []
    for l, lname in enumerate(lnames):
        domain = domainfile[list(domainfile[:,1]).index(lname),-1]
        if domain.split('-')[0] in domaintype:
            keep.append(l)
            dtype.append(domaintype.index(domain.split('-')[0]))
    latentrep = latentrep[keep]
    lnames = lnames[keep]
    bindingkmers = bindingkmers[keep]
else:
    dtype = np.zeros(len(lnames))

difftypes = np.unique(dtype)

latdist = cdist(latentrep, latentrep, 'cosine')
clthresh = float(sys.argv[3])
from sklearn.cluster import AgglomerativeClustering
fig = plt.figure(figsize = (3.5,3.5), dpi = 300)
ax = fig.add_subplot(111)
    


if '--rescalex' in sys.argv:
    add+='_rescalex'
    scaling = []
    for di, ditype in enumerate(difftypes):
        scaling.append(np.array(sys.argv[sys.argv.index('--rescalex')+di+1].split(','),dtype = float))
    


dcolors = ['grey', 'sienna', 'darkolivegreen', 'crimson']
dfitcolors = ['k', 'sienna', 'darkolivegreen', 'crimson']

xlims = []
parameter = []
ylims = []
for di, ditype in enumerate(difftypes):
    
    dmask = dtype == ditype
    print(domaintype[di], 'Currently covered by', np.sum(np.sum(bindingkmers, axis = 0)> 0))
        
    latdistdom = latdist[dmask][:,dmask]
    latbinddom = bindingkmers[dmask]
    
    aggclust = AgglomerativeClustering(n_clusters=None, affinity='precomputed',linkage='complete', distance_threshold=clthresh)
    aggclust.fit(latdistdom)
    seqmotclusters = aggclust.labels_

    clusters = np.unique(seqmotclusters)
    clusterkmers = []
    for clu in clusters:
        clusterkmers.append(np.sum(latbinddom[seqmotclusters == clu], axis = 0)>0)
    clusterkmers = np.array(clusterkmers)
    print(len(clusters))
    x = []
    y = []
    for p in range(500):
        order = np.random.permutation(len(clusters))
        orderk = clusterkmers[order]
        numkmer = np.sum(np.cumsum(orderk, axis = 0) > 0, axis = 1)
        x.append(np.arange(0,len(numkmer)+1))
        y.append(np.append([0],numkmer))

    #if '--domain_split' in sys.argv:
        #y = np.mean(np.array(y), axis = 0)
        #x = np.mean(np.array(x), axis = 0)
    #else:
    x = np.concatenate(x)
    y = np.concatenate(y)
            

    meank = np.mean(np.sum(clusterkmers,axis = 1))
    
    def numkmer(x, a, m):
        y = (1./a)*np.log(m*a*x+1.)
        return y

    def res(par, x, y, mean):
        return numkmer(x,par[0],mean)-y

    def target(x, par, m, y):
        return numkmer(x,par[0],m)-y

    x0 = [meank/4.**7]
    res_soft_l1 = least_squares(res, x0 ,bounds = [[0],[1]], args=(x, y, meank)) #, loss='soft_l1', bounds = [[0,0,0],[np.inf,np.inf, np.inf]],
    params = res_soft_l1.x
    print(domaintype[di],'Params', params[0], meank)

    xlim = least_squares(target, [len(clusters)], args=(params, meank, 4.**7))
    xlim = int(xlim.x[0])
    print( 'Xlim', xlim)
    
    parameter.append([params[0], meank])
    
    
    if '--rescalex' in sys.argv:
        scaling[di][0] += len(clusters)
        if xlim < scaling[di][0]:
            maxc = xlim
        else:
            maxc = scaling[di][0]*100./scaling[di][1]
        xtest = np.arange(1,maxc)
        xlims.append([np.around((float(len(clusters))/scaling[di][0])*scaling[di][1],1), np.around((float(maxc)/scaling[di][0])*scaling[di][1],1)])
    else:
        maxc = xlim
        xtest = np.append(np.concatenate([np.arange(1,10)*10**i for i in range(int(np.log10(xlim)))]),[xlim])
        xlims.append([len(clusters), xlim])
        
    ylims.append(numkmer(np.array([1,maxc]), params[0], meank))
    ytest = numkmer(xtest, params[0], meank)

    # probability of overlap
    def poverlap(lap):
        return 1.-np.exp(-lap)
    
    a = np.around(poverlap(params[0])*100.,4)

    if '--domain_split' in sys.argv:
        yn = []
        for sx in np.unique(x):
            yn.append(np.mean(y[x==sx]))
        y = np.array(yn)
        x = np.unique(x)
        

    if '--rescalex' in sys.argv:
        xtest = (xtest/scaling[di][0])*scaling[di][1]
        x = (x/scaling[di][0])*scaling[di][1]
        
    ax.fill_between(xtest,0,ytest, color = dcolors[di], alpha = 0.2)#, label = str(int(xlim-len(clusters)))+' additional cluster')
    ax.scatter(x,y, color = dcolors[di], alpha = 0.2, label = domaintype[di]+' measured clusters')
    ax.plot(xtest, ytest, ls = '--', c = dfitcolors[di], label = r'$mean={0}, p_\cap={1}\%$'.format(str(np.around(meank,1)), a)) # r'$k_{'+domaintype[di]+'}=1/a \cdot ln(a \cdot m \cdot x+1)$'+'\n'+
    #ax.plot([0, xlims[-1][0]], [y[-1],y[-1]], c = dcolors[di], ls = ':')
    if '--rescalex' in sys.argv:
        ax.text(xlims[-1][0]+5, y[-1], '('+str(xlims[-1][0])+'%,'+str(int(y[-1]*100./4.**7))+'%)')
    else:
        ax.text(xlims[-1][0]+5,y[-1], '('+str(xlims[-1][0])+','+str(int(y[-1]*100./4.**7))+'%)')
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_ylabel('7-mers bound by RBP')

if '--rescalex' in sys.argv:
    ax.set_xlabel('Measurable sequence-specificity clusters')
    ax.set_xticks(np.linspace(0,100,6))
    ax.set_xticklabels(["0%","20%","40%","60%","80%",'100%'])
else:
    ax.set_xlabel('Sequence-specificity clusters')

if '--logx' in sys.argv:
    ax.set_xscale('log')
else:
    ax.set_xlim([0,np.amax(xlims)])

if '--logy' in sys.argv:
    ax.set_yscale('log')
    ax.set_yticks([16.384,163.84,1638.4, 4**7])
    ax.set_ylim([16,4**7+1])
    ax.set_yticklabels(['0.1%','1%', '10%', '100%'])
else:
    ax.set_yticks(np.linspace(0,1.,6)*4**7)
    ax.set_yticklabels(['0%','20%','40%','60%', '80%', '100%'])
    ax.set_ylim([0,4**7])

for di, xlim in enumerate(xlims):
    print( ylims[di][-1]/4.**7)
    if ylims[di][-1]/4.**7 < 0.99:
        perach = np.around((ylims[di][-1]*100.)/4.**7,1)
        ax.text(xlim[-1], ylims[di][-1], str(perach)+'%')
    else:
        if '--rescalex' in sys.argv:
            ax.text(xlim[-1], ylims[di][-1], str(xlim[-1])+'%')
        else:
            ax.text(xlim[-1], ylims[di][-1], str(xlim[-1]))
            
        
        
ax.grid()
#ax.set_xticklabels(np.concatenate(xlims))



ax.legend(loc = 0, prop={'size':6})

if '--savefig' in sys.argv:
    fig.savefig('K-mer_coverage_seqmotclusters'+add+'-cut'+str(clthresh)+'.jpg', dpi = 300, bbox_inches = 'tight')

plt.show()










