import numpy as np
import sys, os
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
# cluster on latent distances
from scipy.optimize import least_squares
from scipy import stats
from statsmodels.stats.multitest import multipletests


zscores = np.genfromtxt(sys.argv[1], dtype = float)[:,1:]
znames = open(sys.argv[1],'r').readline().strip().split()[1:]

bindingkmers = []
kmerpvalues = []
for z, zname in enumerate(znames):
    pval = 2.*(1 - stats.norm.cdf(zscores[:,z]))
    p_valbind = multipletests(pval, alpha=0.001, method='fdr_bh', is_sorted=False, returnsorted=False)
    kmerpvalues.append(p_valbind[1])
    bindingkmers.append(p_valbind[0])
bindingkmers = np.array(bindingkmers)
kmerpvalues =np.array(kmerpvalues)

domains = ['RRM', 'KH']
domainfile = np.genfromtxt(sys.argv[2], dtype = str)
doms = []
for z in range(len(domainfile)):
    domainfile[z][-1] = '-'.join(np.unique(domainfile[z][-1].split('-')))



x = np.array([0.01,0.001])
y = [np.sum(np.sum(kmerpvalues<=p, axis = 0) > 0) for p in x]
yrrm = [np.sum(np.sum(kmerpvalues[np.isin(znames,domainfile[domainfile[:,-1]=='RRM',1])]<=p, axis = 0) > 0) for p in x]
ykh = [np.sum(np.sum(kmerpvalues[np.isin(znames,domainfile[domainfile[:,-1]=='KH',1])]<=p, axis = 0) > 0) for p in x]
drrm= np.sum(kmerpvalues[np.isin(znames,domainfile[domainfile[:,-1]=='RRM',1])]<=0.01, axis = 1) 
dkh= np.sum(kmerpvalues[np.isin(znames,domainfile[domainfile[:,-1]=='KH',1])]<=0.01, axis = 1)
y = 100*np.array(y)/(4**7)
yrrm = 100*np.array(yrrm)/(4**7)
ykh = 100*np.array(ykh)/(4**7)

fig = plt.figure(figsize=(3,3), dpi = 300)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.bar(np.arange(len(y))-0.45,y, width=0.3, color = 'dimgrey', alpha = 0.8, ls='-', lw = 1, edgecolor = 'dimgrey', label = 'All RBPs')
ax.bar(np.arange(len(y))-0.15,yrrm,width=0.3, color = 'purple', alpha = 0.4, ls='-', lw =1, edgecolor = 'purple', label = 'RRM RBPs')
ax.bar(np.arange(len(y))+0.15,ykh,width=0.3, color = 'firebrick', alpha = 0.8, ls='-', lw = 1, edgecolor = 'firebrick', label='KH RBPs')
ax.set_xticks(np.arange(len(y)))
ax.set_xticklabels(x)
ax.set_ylim([0,100])
ax.set_xlabel('FDR(Zscore)')
ax.set_ylabel('% recognized 7-mers')
ax.legend(fontsize = 8)
fig.savefig(sys.argv[3], bbox_inches = 'tight', dpi = 500)

figh = plt.figure(figsize=(3,3), dpi = 300)
bins = np.arange(0,380,20)
axh = figh.add_subplot(111)
axh.spines['top'].set_visible(False)
axh.spines['right'].set_visible(False)
axh.hist(drrm, color = 'purple', bins = bins,label='RRM RBPs', alpha = 0.4, histtype = 'stepfilled')
axh.hist(dkh, color = 'firebrick',bins = bins,label='KH RBPs', alpha = 0.8, histtype = 'stepfilled')
axh.hist(drrm, color = 'purple', bins = bins, histtype = 'step')
axh.hist(dkh, color = 'firebrick',bins = bins, histtype = 'step')
axh.legend(fontsize = 8)
axh.set_xlabel('# 7-mers with FDR<0.01')
axh.set_ylabel('# RBPs')
figh.savefig(os.path.splitext(sys.argv[3])[0]+'_indvnumbers'+os.path.splitext(sys.argv[3])[1], bbox_inches = 'tight', dpi = 500)



