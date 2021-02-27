import numpy as np
import sys, os
from scipy.stats import ttest_ind, ttest_1samp, norm
from statsmodels.stats.multitest import multipletests


sets = sys.argv[1].split(':')
nsets = int(sys.argv[2])
# load original file
ofile = np.genfromtxt(sys.argv[3], dtype = str)
pcut = float(sys.argv[5])
mcut = float(sys.argv[4])



zscores = ofile[:,1].astype(float)
rbps = ofile[:,0]
pvalues = ofile[:,2].astype(float)
# load all shuffled files
szscores = []
for s in range(nsets+1):
    if os.path.isfile(sets[0]+str(s)+sets[1]):
        sfile = np.genfromtxt(sets[0]+str(s)+sets[1], dtype = str)
        szscore = sfile[:,1].astype(float)
        sbps = sfile[:,0]
        if np.array_equal(sbps, rbps):
            szscores.append(szscore)
        else:
            print sbps, rbps
    else:
        print sets[0]+str(s)+sets[1], 'missing'
        

def Ztest(X1, me):
    sd = np.std(X1)
    m1 = np.mean(X1)
    z = (me - m1)/sd
    pval = 2.*(1 - norm.cdf(abs(z)))
    return round(z, 3), pval


szscores = np.array(szscores)
zofz = []
pofp = []
for r, rbp in enumerate(rbps):
    z, pv = Ztest(szscores[:,r], zscores[r])
    zofz.append(z)
    pofp.append(pv)
pofp = np.array(pofp)
zofz = np.array(zofz)

multipfdr = multipletests(pvalues, alpha = mcut, method = 'fdr_bh')
#multipofpfdr = multipletests(pofp, alpha = pcut, method = 'fdr_bh')

mpval_adaptfdr = multipfdr[1]
mpofp_adaptfdr = pofp
maskfdr = multipfdr[0]*(pofp<=pcut)*(zscores*zofz > 0)



outobj = open(os.path.splitext(sys.argv[3])[0]+'-fdr'+str(pcut)+'_zofz'+str(mcut)+'.pvals', 'w')
outobj.write('#Z ZofZ pval Pofz fdrp fdrpofz fdraccept\n')
for r, rbp in enumerate(rbps):
    outobj.write(rbp+' '+str(zscores[r])+' '+str(zofz[r])+' '+str(mpval_adaptfdr[r])+' '+str(mpofp_adaptfdr[r])+' '+str(maskfdr[r])+'\n')






