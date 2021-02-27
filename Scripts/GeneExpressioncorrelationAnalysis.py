import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from scipy.stats import ttest_ind, ttest_1samp, ttest_ind_from_stats
from scipy.spatial.distance import cdist
from scipy import stats
import time
from scipy.optimize import Bounds
from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest

def check(var):
    if var == 'true' or var == 'True':
        return True
    else:
        return False

def twoSampZ(X1, X2):
    sd1,sd2 = np.std(X1),np.std(X2)
    n1,n2 = float(len(X1)), float(len(X2))
    m1,m2 = np.mean(X1),np.mean(X2)
    pooledSE = np.sqrt(sd1**2/n1 + sd2**2/n2)
    z = (m1 - m2)/pooledSE
    pval = 2.*(1 - stats.norm.cdf(abs(z)))
    return round(z, 3), pval

def MWUtest(X1, X2):
    ust, pv = stats.mannwhitneyu(X1, X2, alternative= 'two-sided')
    pv = pv/2.
    n1,n2 = float(len(X1)), float(len(X2))
    z = (ust - (n1*n2/2.))/np.sqrt(n1*n2*(n1+n2+1.)/12.)
    return round(z,3), pv

def correlationanalysis(bindingmat, protexpression, expressions, consider = 'all', method = 'correlation', conditions = None, cutbind = 1., cutmethod = 'cutoff', combined = False):
    
    if consider == 'expressed':
        median = np.median(protexpression, axis = 1)
        consider = (protexpression.T - median).T > 0
    else:
        consider = np.ones(np.shape(protexpression)) == 1
    
    nprots = len(protexpression)
    if combined:
        protexpression = (protexpression.T/(np.amax(protexpression, axis =1)-np.amin(protexpression, axis =1))).T
        dprotn = (nprots*(nprots-1))/2
        dprotexpression = []
        dbindingmat = []
        for p, protexpress in enumerate(protexpression):
            for q in range(p+ 1, len(protexpression)):
                qrotexpress = protexpression[q]
                dprotexpression.append(protexpress*qrotexpress)
                bmat = -np.ones(len(bindingmat[p]))
                bmat[(bindingmat[p]>=cutbind)*(bindingmat[q]>=cutbind)] = cutbind
                dbindingmat.append(bmat)
        bindingmat = np.array(dbindingmat)
        protexpression = np.array(dprotexpression)
    
    
    if conditions is None:
        correlationmat = []
        for p, protexpress in enumerate(protexpression):
            correlationmat.append(1. - cdist([protexpress[consider[p]]], expressions[:, consider[p]], 'correlation')[0])
        correlationmat = np.nan_to_num(np.array(correlationmat))
    else:
        condchange = np.zeros(len(expression))
        for g, gene in enumerate(expression):
            if np.std(gene) != 0:
                z, pv = MWUtest(gene[conditions == 1], gene[conditions == 0])
                condchange[g] = z
        rbpcondchange = np.zeros(len(protexpression))
        for g, gene in enumerate(protexpression):
            if np.std(gene) != 0:
                z, pv = MWUtest(gene[conditions == 1], gene[conditions == 0])
                if pv < 0.075:
                    rbpcondchange[g] = 1.
                    
        correlationmat = np.dot(rbpcondchange.reshape(-1,1),condchange.reshape(1,-1))
        print np.shape(correlationmat), np.shape(bindingmat)
    
    if cutmethod == 'fdr':
        for b, bmat in enumerate(bindingmat):
            floc, explam = stats.norm.fit(bmat)
            p_valbind = stats.norm.sf(bmat, loc = floc, scale = explam)
            p_valbind = multipletests(p_valbind, alpha=cutbind, method='fdr_bh', is_sorted=False, returnsorted=False)
            bindingmat[b] = -np.log10(p_valbind[1])
        cutbind = -np.log10(cutbind)
    elif cutmethod == 'best':
        for b, bdmat in enumerate(bindingmat):
            bindingmat[b] = bdmat/np.sort(np.unique(bdmat))[-int(cutbind)]
        cutbind = 1.
    elif cutmethod == 'top':
        for b, bdmat in enumerate(bindingmat):
            bindingmat[b] = bdmat/np.sort(bdmat)[-int(cutbind)]
        cutbind = 1.
    elif cutmethod == 'optimal':
        for b, bdmat in enumerate(bindingmat):
            cormat = correlationmat[b]
            diffset = [0]
            unvalues = np.sort(np.unique(bdmat))
            if np.std(cormat) > 0:
                for i in range(0,min(len(np.unique(bdmat))-1,int(cutbind))):
                    bmattest = unvalues[-i-1]
                    if len(cormat[bdmat >= bmattest]) > 0.1*len(cormat):
                        break
                    if len(cormat[bdmat >= bmattest]) > 10:
                        diff, pd = MWUtest(cormat[bdmat >= bmattest], cormat[bdmat < bmattest])
                    else:
                        diff = 0
                    diffset.append(np.absolute(diff))
                ldif = np.argmax(diffset)
            else:
                ldif = 1
            bindingmat[b] = bdmat/unvalues[-ldif]
        cutbind = 1.
    
    
    if combined:
        dconsider = []
        for p, cons in enumerate(consider):
            for q in range(p+ 1, len(consider)):
                dconsider.append(cons*consider[q])
        consider = np.array(dconsider)
    

    
    z_scores = []
    p_values = []
    mean_corr = []
    max_corr = []
    for p, pstats in enumerate(correlationmat):
        if len(correlationmat[p, bindingmat[p] >= cutbind])> 0 and len(correlationmat[p, bindingmat[p] < cutbind])> 0:
            if method == 'correlation':
                if np.std(correlationmat[p]) != 0:
                    z, pv = MWUtest(correlationmat[p,bindingmat[p] >= cutbind], correlationmat[p, bindingmat[p] < cutbind])
                else:
                    z,pv = 0.,1.
            elif method == 'binding':
                corrmax = np.amax(np.absolute(pstats))
                corrcut = np.amax(np.absolute(pstats))
                zs = []
                ps = []
                corcuts = []
                bound = []
                for i in range(49):
                    corrcut -= 0.02*corrmax
                    if np.sum(pstats >= corrcut) > 0:
                        A, a, B, b = np.sum(pstats >= corrcut), np.sum(bindingmat[p][pstats >= corrcut]>=cutbind), np.sum(pstats < corrcut), np.sum(bindingmat[p][pstats < corrcut]>=cutbind) 
                        z, pv = stats.fisher_exact([[a,b],[A-a, B-b]], alternative = 'greater')
                        zs.append(z)
                        ps.append(pv)
                        corcuts.append(corrcut)
                        bound.append(np.around([a/float(A),b/float(B)],2))
                    if np.sum(pstats <= -corrcut) > 0:
                        A, a, B, b = np.sum(pstats <= -corrcut),np.sum(bindingmat[p][pstats <= -corrcut]>=cutbind), np.sum(pstats > -corrcut), np.sum(bindingmat[p][pstats > -corrcut]>=cutbind)
                        z, pv = stats.fisher_exact([[a,b],[A-a, B-b]], alternative = 'greater')
                        zs.append(z)
                        ps.append(pv)
                        corcuts.append(-corrcut)
                        bound.append(np.around([a/float(A),b/float(B)],2))
                minp = np.argmin(ps)
                pv, z, corrcut = ps[minp], zs[minp], corcuts[minp]
                nbmat = np.zeros(len(bindingmat[p]))
                nbmat[(np.sign(corrcut)*pstats >= np.sign(corrcut)*corrcut) * (bindingmat[p] >= cutbind)] = cutbind
                bindingmat[p]= nbmat
            else:
                print 'Analysis method not understood'
                sys.exit()
            
            means = np.around([np.mean(correlationmat[p, bindingmat[p] >= cutbind]), np.mean(correlationmat[p, bindingmat[p] < cutbind])],2)
            maxc = np.around([np.percentile(correlationmat[p, bindingmat[p] >= cutbind], [5,50,95]), np.percentile(correlationmat[p, bindingmat[p] < cutbind], [5,50,95])],2)
        else:
            z, pv = 0.,1.
            means = [0.,0.]
            maxc = [[0,0,0],[0,0,0]]
            
        max_corr.append(maxc)
        mean_corr.append(means)
        z_scores.append(z)
        p_values.append(pv)
    
    z_scores = np.nan_to_num(z_scores)
    p_values = np.array(p_values)
    p_values[np.isnan(p_values)] = 1.
    explaindata = np.amax(z_scores)
    return explaindata, correlationmat, z_scores, p_values, np.array(mean_corr), np.array(max_corr), bindingmat>=cutbind
    


    




if __name__ == '__main__':
    expressionfile = np.load(sys.argv[1])
    
    expression = expressionfile['pkm']
    tissues = expressionfile['tissues']
    genes = expressionfile['genes']
    #for g, gene in enumerate(genes):
        #genes[g] = gene.split('.')[0]
    
    
    bindingfile = np.load(sys.argv[2])
    bindingmat = bindingfile['bindingmat']
    rbpnames = bindingfile['rbpnames']
    genenames = bindingfile['genenames']
    
    
    if '--outname' in sys.argv:
        outname = sys.argv[sys.argv.index('--outname')+1]
    else:
        outname = 'out'
    
    
    if not np.array_equal(genes, genenames):
        print "align genes and genenames"
        print '...'
        uniquegenes, genemask = np.unique(genes, return_index = True)
        uniquegenenames, genenamemask = np.unique(genenames, return_index = True)
        
        genemask = genemask[np.isin(uniquegenes, uniquegenenames)]
        genenamemask = genenamemask[np.isin(uniquegenenames, uniquegenes)]
        
        expression = expression[genemask]
        genes = genes[genemask]
        bindingmat = bindingmat[:, genenamemask]
        genenames = genenames[genenamemask]
        print np.shape(expression), np.shape(bindingmat)
        
        if not np.array_equal(genes, genenames):
            geneorder = []
            genenameorder = []
            for g, gene in enumerate(genenames):
                gn = gene.split()[1].split(':')[1].split('.')[0]
                if gn in genes:
                    geneorder.append(list(genes).index(gn))
                    genenameorder.append(g)    
        
            expression = expression[geneorder]
            bindingmat = bindingmat[:, genenameorder]
        
        print 'done', np.shape(expression), np.shape(bindingmat)
    
    rbpfile = np.load(sys.argv[3])
    rbpexpression = rbpfile['pkm']
    stissues = rbpfile['tissues']
    rbps = rbpfile['genes']
    
    
    if not np.array_equal(rbpnames, rbps):
        print 'rbps and rbpnames need to be aligned'
        rorder = []
        border = []
        for r, rbp in enumerate(rbpnames):
            if rbp in rbps:
                border.append(r)
                rorder.append(list(rbps).index(rbp))
        rbps = rbps[rorder]
        rbpexpression = rbpexpression[rorder]
        rbpnames = rbpnames[border]
        bindingmat = bindingmat[border]
        print np.shape(bindingmat), np.shape(rbpnames)
    
    if not np.array_equal(stissues, tissues):
        print 'align tissue'
        sort = []
        for t, tis in enumerate(tissues):
            sort.append(list(stissues).index(tis))
        stissues = stissues[sort]
        rbpexpression = rbpexpression[:, sort]
        print np.shape(rbpexpression)
        
    
    if '--tissuemean' in sys.argv:
        outname += '-tissuemean'
        tfile = open(sys.argv[sys.argv.index('--tissuemean')+1], 'r').readlines()
        tissueid = []
        for l, line in enumerate(tfile):
            if l> 0:
                line = line.split('\t')
                tissueid.append([line[8], line[10]])
        tissueid = np.array(tissueid)
        uniquetissue, tn = np.unique(tissueid[:, 1], return_index = True)
        nrbpexpression = np.zeros((len(rbps), len(uniquetissue)))
        nexpression = np.zeros((len(genes), len(uniquetissue)))
        
        locations = []
        for u, utiss in enumerate(uniquetissue):
            tissueids = tissueid[tissueid[:,1] == utiss, 0]
            loc = np.where(np.isin(tissues, tissueids))[0]
            locations.append(loc)
            nrbpexpression[:, u] = np.mean(rbpexpression[:, loc], axis = 1)
            nexpression[:, u] = np.mean(expression[:, loc], axis = 1)
        locations = np.array(locations).T
        
        # filter genes with bad correlation between replicates R> 0.7
        expcosine = np.sum(expression[:, locations[0]].T/np.sqrt(np.sum(expression[:, locations[0]]**2, axis =1)) * expression[:, locations[1]].T/np.sqrt(np.sum(expression[:, locations[1]]**2, axis =1)), axis = 0)
        print len(expcosine), np.amax(expcosine), np.amin(expcosine)
        expmask = expcosine > 0.6
        
        rbpcosine = np.sum(rbpexpression[:, locations[0]].T/np.sqrt(np.sum(rbpexpression[:, locations[0]]**2, axis =1)) * rbpexpression[:, locations[1]].T/np.sqrt(np.sum(rbpexpression[:, locations[1]]**2, axis =1)), axis = 0)
        print len(rbpcosine), np.amin(rbpcosine)
        rbpexpmask = rbpcosine > 0.7
        
        
        
        expression = nexpression[expmask]
        bindingmat = bindingmat[rbpexpmask][:, expmask]
        rbpexpression = nrbpexpression[rbpexpmask]
        tissues = stissues = uniquetissue
        
        genes = genenames = genes[expmask]
        rbps = rbps[rbpexpmask]
        rbpnames = rbpnames[rbpexpmask]
        
        print 'mean tissue', np.shape(expression), np.shape(rbpexpression)

    if '--lognorm_expression' in sys.argv:
        for r, exp in enumerate(expression):
            expression[r][exp <= 0] = np.amin(exp[exp > 0])
        expression = np.log2(expression/np.mean(expression, axis = 1)[:,None])

    if '--lognorm_rbpexpression' in sys.argv:
        for r, rbpexp in enumerate(rbpexpression):
            rbpexpression[r][rbpexp <= 0] = np.amin(rbpexp[rbpexp > 0])
        rbpexpression = np.log2(rbpexpression/np.mean(rbpexpression,axis = 1)[:, None])
        


    if '--select_tissue' in sys.argv:
        
        tissuefile = np.genfromtxt(sys.argv[sys.argv.index('--select_tissue')+1], dtype = str, delimiter = '\t')
        outname += '-on'+os.path.splitext(os.path.split(sys.argv[sys.argv.index('--select_tissue')+1])[1])[0]
        medianexpression = (np.amax(expression, axis = 1)+np.amin(expression, axis = 1))/2.
        medianrbpexpression = (np.amax(rbpexpression, axis = 1)+np.amin(rbpexpression, axis = 1))/2.
        tmask = np.isin(tissues, tissuefile)
        tissues = tissues[tmask]
        expression = expression[:, tmask]
        rbpexpression = rbpexpression[:, tmask]
        
        expmask = np.amax(expression - medianexpression[:, None], axis = 1) > 0
        rbpexpmask = np.amax(rbpexpression - medianrbpexpression[:, None], axis = 1) > 0
        
        expression = expression[expmask]
        bindingmat = bindingmat[rbpexpmask][:, expmask]
        rbpexpression = rbpexpression[rbpexpmask]
        genes = genenames = genes[expmask]
        rbps = rbps[rbpexpmask]
        rbpnames = rbpnames[rbpexpmask]
        print 'selected tissue', np.shape(expression), np.shape(rbpexpression)

    if '--two_conditions' in sys.argv:
        conditions = []
        for t in tissues:
            if 'A' in t:
                conditions.append(1)
            else:
                conditions.append(0)
        conditions = np.array(conditions, dtype = int)
    else:
        conditions = None
    
    
    if '--correlationanalysis' in sys.argv:
        corrmethod = sys.argv[sys.argv.index('--correlationanalysis')+1] # expressed or all
        bindcut = float(sys.argv[sys.argv.index('--correlationanalysis')+2])
        double = check(sys.argv[sys.argv.index('--correlationanalysis')+3])
        cutmethod = sys.argv[sys.argv.index('--correlationanalysis')+4] #fdr, cut, top
        method = sys.argv[sys.argv.index('--correlationanalysis')+5] #correlation, binding
        
        explaindata, correlationmat, zscores, p_values, mediancorr,maxcorr, bindingmat = correlationanalysis(bindingmat, rbpexpression, expression, consider = corrmethod, method = method, conditions = conditions, cutbind = bindcut, cutmethod = cutmethod, combined = double)
        
        
        if double:
            nrbpnames = []
            for p, prot in enumerate(rbpnames):
                for q in range(p+1,len(rbpnames)):
                    nrbpnames.append(prot+'--'+rbpnames[q])
            rbpnames = np.array(nrbpnames)
            
        multipdfr = multipletests(p_values, alpha = 0.1, method = 'fdr_bh')       
        
        mtpcut = 0.01
        
        if len(mediancorr[multipdfr[1] < mtpcut]) > 0:
            print len(mediancorr[multipdfr[1] < mtpcut])
            for g in np.argsort(-multipdfr[1])[np.where(-np.sort(-multipdfr[1]) < mtpcut)[0]]:
                print rbpnames[g], g, zscores[g], p_values[g], multipdfr[1][g], -np.diff(mediancorr[g]), maxcorr[g][0], maxcorr[g][1], np.sum(bindingmat[g])
        
        print 'MinP', np.amin(p_values)*float(len(p_values))
        
        if '--npz' in sys.argv:
            np.savez_compressed(os.path.splitext(sys.argv[2])[0]+'-'+outname+'_'+method+'_'+corrmethod+'_'+cutmethod+str(bindcut)+str(double)+'_correlationanalysis.npz', correlationmat = correlationmat, bindingmat = bindingmat, rbpnames=rbpnames, genes=genes)
        
        np.savetxt(os.path.splitext(sys.argv[2])[0]+'-'+outname+'_'+method+'_'+corrmethod+'_'+cutmethod+str(bindcut)+str(double)+'_correlationanalysis.txt', np.concatenate([[rbpnames], [zscores], [p_values], -np.around(np.diff(mediancorr, axis = 1),2).T, [np.sum(bindingmat, axis = 1).astype(int)],[multipdfr[1]],[multipdfr[1]< mtpcut]], axis = 0).T.astype(str), fmt = '%s')
        
        
    
    
    
    
