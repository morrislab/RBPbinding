import numpy as np
import sys, os
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import glob

if '--determine_next' in sys.argv:

    species = np.sort(glob.glob('*/')) #genfromtxt(sys.argv[1], dtype = str)
    if '--take_sample' in sys.argv:
        sampsize = int(sys.argv[sys.argv.index('--take_sample')+1])
        species = species[np.random.permutation(len(species))[:sampsize]]
    latentrep = sys.argv[1]
    domainfile = np.genfromtxt(sys.argv[2], dtype = str, delimiter = ' ', comments = None)
    latentmeas = np.load(sys.argv[3])
    measprots = latentmeas['names']
    measlatent = latentmeas['profiles']
    print np.shape(measlatent)
    
    if '--confidence' in sys.argv:
        confcutset = [float(sys.argv[sys.argv.index('--confidence')+1])]
    else:
        confcutset = [0.127,0.2, 0.33]

    
    latent = []
    latentnames = []
    for s, spec in enumerate(species):
        slat = np.load(spec+'/'+latentrep)
        latent.append(slat['profiles'])
        latentnames.append(slat['names'])

    latentnames = np.concatenate(latentnames)
    latent = np.concatenate(latent, axis = 0)
    print len(latent)
    
    for d, dom in enumerate(domainfile):
        domainfile[d,-1] = dom[-1].split('-')[0]
    
    lsort = np.argsort(latentnames)
    latentnames = latentnames[lsort]
    latent = latent[lsort]
    dsort = np.argsort(domainfile[:,0])
    domainfile = domainfile[dsort]
    if not np.array_equal(domainfile[:,0], latentnames):
        print 'domains dont match'
        lsort = np.isin(latentnames, domainfile[:,0])
        latentnames = latentnames[lsort]
        latent = latent[lsort]
        dsort = np.isin(domainfile[:,0], latentnames)
        domainfile = domainfile[dsort]
        if not np.array_equal(domainfile[:,0], latentnames):
            print len(domainfile[:,0]), len(latentnames), len(np.unique(latentnames)), len(np.unique(domainfile[:,0])), domainfile[0,0], latentnames[0]
            sys.exit()
    
    
    allnames = np.copy(latentnames)
    fullen = float(len(latentnames))
    for co, confcut in enumerate(confcutset):
        latentnames = np.copy(allnames)
        cluster = -np.ones(int(fullen), dtype = int)
        
        # go trough measured proteins
        print np.shape(latent),np.shape(measlatent)
        latdist = cdist(measlatent, latent, 'cosine')
        mlatdist = cdist(measlatent,measlatent, 'cosine')
        ci = 0
        while True:
            cnum = np.sum(latdist<=confcut, axis = 1)
            argmax = np.argmax(cnum)
            
            if cnum[argmax] == 0:
                break
            else:
                mprots = mlatdist[argmax]<=confcut
                
                cprots = latdist[argmax]<=confcut
                cluster[np.isin(allnames,latentnames[cprots])] = ci
                print ci, np.sum(cluster == ci), np.sum(mprots)
                
                mlatdist = mlatdist[~mprots][:,~mprots]
                latdist = latdist[~mprots][:,~cprots]
                latentnames = latentnames[~cprots]
                ci += 1
                if len(mlatdist) == 0:
                    break
        measclusts = np.copy(ci)
        
        leftmask = np.isin(allnames,latentnames)
        latdist = cdist(latent[leftmask], latent[leftmask], 'cosine')
        while True:
            cnum = np.sum(latdist<=confcut,axis = 1)
            argmax = np.argmax(cnum)
            cprots = latdist[argmax]<=confcut
            
            print ci, np.sum(cprots)
            
            cluster[np.isin(allnames,latentnames[cprots])] = ci
            ci += 1
            latdist = latdist[~cprots][:,~cprots]
            latentnames = latentnames[~cprots]
            
            if cnum[argmax] == 1:
                print len(latentnames), 'clusters with only 1'
                lcluster = np.arange(ci, ci +len(latentnames), dtype = int)
                cluster[np.isin(allnames,latentnames)] = lcluster
                break
        
        measured = (cluster < measclusts).astype(int)
        np.savetxt('JPLEcosineclusters_cut'+str(confcut)+'.txt', np.array([allnames,cluster,domainfile[:,-1], measured]).T, fmt = '%s')

else:
    recoveredlist = []
    if ',' in sys.argv[1]:
        files = sys.argv[1].split(',')
        cuts = sys.argv[2].split(',')
    else:
        files = [sys.argv[1]]
        cuts = [sys.argv[2]]
    
    for cfile in files:
        print cfile
        clusterfile = np.genfromtxt(cfile, dtype = str, delimiter = ' ')
        clusterfile = clusterfile[((clusterfile[:,-2] == 'KH') | (clusterfile[:,-2] == 'RRM'))]
        clusters, cindex, recovered = np.unique(clusterfile[:,1].astype(int), return_index = True, return_counts = True)
        rrmfraction = np.zeros(len(clusters))
        clusterdomains = clusterfile[:, 2][cindex]
        rrmfraction[(recovered == 1) & (clusterdomains== 'RRM')] = 1
        for c, cl in enumerate(clusters):
            if c%1000 ==0:
                print c, len(clusters)
                
            if recovered[c] > 1:
                clmask = clusterfile[:,1].astype(int)==cl
                rrmfraction[c] = np.sum(clusterfile[clmask,2] == 'RRM')/float(recovered[c])
        
        recoveredlist.append([rrmfraction, clusterfile[cindex,-1]])
        

    fig = plt.figure(figsize = (4,3.*len(recoveredlist)))

    for s, csets in enumerate(recoveredlist):
        fracmix = np.sum((csets[0] != 0) & (csets[0] != 1.))/float(len(csets[0]))
        print 'False unification', fracmix
        
        fracmixm = np.sum((csets[0][csets[1] == '1'] != 0) & (csets[0][csets[1] == '1'] != 1.))/float(len(csets[0][csets[1] == '1']))
        print 'False unification measured', fracmixm
        
        ax = fig.add_subplot(len(recoveredlist),1,s+1)
        ax.hist(csets[0], bins = 20, color = 'grey', alpha = 0.5, label = 'Unmeasured: '+str(round(fracmix,3)) )
        ax.hist(csets[0][csets[1] == '1'], bins = 20, alpha = 0.5, color = 'steelblue', label = 'Measured: '+str(round(fracmixm,3)))
        

        ax.set_yscale('symlog')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticks([1,10,100,1000,10000])
        ax.set_yticklabels([1,10,100,1000,10000])
        ax.legend()
        ax.set_ylabel('Number clusters')
        ax.grid(axis = 'y')
        ax.set_title('FP rate cut: '+ str(cuts[s]))
    ax.set_xlabel('Fraction RRMs')

    if '--savefig' in sys.argv:
        outname = 'JPLEclustering_RRMdistribution'+'-'+'-'.join(np.array(cuts).astype(str))+'.jpg'
        outname = outname.replace(' ', '-')
        print outname
        fig.savefig(outname, bbox_inches = 'tight', dpi = 300)
    
    plt.show()

