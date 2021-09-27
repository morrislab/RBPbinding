import numpy as np
import sys, os
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import glob

# split KH, RRM

if '--determine_next' in sys.argv:

    species = np.sort(glob.glob('*/')) #genfromtxt(sys.argv[1], dtype = str)
    if '--take_sample' in sys.argv:
        sampsize = int(sys.argv[sys.argv.index('--take_sample')+1])
        species = species[np.random.randint(0,len(species), sampsize)]
    latentrep = sys.argv[1]
    
    
    if '--confidence' in sys.argv:
        confcutset = [float(sys.argv[sys.argv.index('--confidence')+1])]
    else:
        confcutset = [0.127,0.2,0.3]

    if '--maxprots' in sys.argv:
        maxnum = int(sys.argv[sys.argv.index('--maxprots')+1])
    else:
        maxnum = 200

    latent = []
    latentnames = []
    for s, spec in enumerate(species):
        slat = np.load(spec+'/'+latentrep)
        latent.append(slat['profiles'])
        latentnames.append(slat['names'])

    latentnames = np.concatenate(latentnames)
    latent = np.concatenate(latent, axis = 0)
    print len(latent)
    
    add = ''
    if '--domain_only' in sys.argv:
        domainfile = np.genfromtxt(sys.argv[sys.argv.index('--domain_only')+1], dtype = str, delimiter = ' ', comments = None)
        domaintype = sys.argv[sys.argv.index('--domain_only')+2]
        add += '-domain'+domaintype
        keep = []
        for d, pd in enumerate(domainfile):
            if domaintype in pd[-1]:
                keep.append(pd[0])
        dmask = np.isin(latentnames, keep)
        latentnames = latentnames[dmask]
        latent = latent[dmask]
        print domaintype, len(latent)
    
    allnames = np.copy(latentnames)
    fullen = float(len(latentnames))
    cluster = -np.ones(int(fullen), dtype = int)


    c = 0
    if '--removeconfident' in sys.argv:
        latentclose = np.genfromtxt(sys.argv[sys.argv.index('--removeconfident')+1], dtype = str)
        
        confident = latentclose[latentclose[:,-1].astype(float)<confcutset[0],0]
        print np.sum(np.isin(latentnames, confident))
    
        latentmask = ~np.isin(latentnames, confident)
    
        cluster[~latentmask] = 0
        latentnames = latentnames[latentmask]
        latent = latent[latentmask]
        print len(latentnames)
        recovered = [np.sum(~latentmask)/fullen]
        c += 1
    
    

        
    
    
    if '--only_within' in sys.argv:
        latentstats = np.genfromtxt(sys.argv[sys.argv.index('--only_within')+1], dtype = str)
        withincut = float(sys.argv[sys.argv.index('--only_within')+2])
        add += '-within'+str(withincut)
        confident = latentstats[latentstats[:,-1].astype(float)<withincut,0]
    else:
        confident = np.copy(latentnames)
        
        
    
    
    recoveredlist = []
    
    latdistin = cdist(latent, latent, 'cosine')
    latentnamesin = np.copy(latentnames)
    for i, confcut in enumerate(confcutset):
        latdist = np.copy(latdistin)
        latentnames = latentnamesin
        confiset = np.copy(confident)
        if '--removeconfident' in sys.argv and i != 0:
            close = latentclose[latentclose[:,-1].astype(float)<confcut,0]
            print len(latentnames), np.sum(np.isin(latentnames, close))
        
            latentmask = ~np.isin(latentnames, close)
        
            cluster[np.isin(allnames, close)] = 0
            latentnames = latentnames[latentmask]
            latdist = latdist[latentmask][:,latentmask]
            print len(latentnames)
            recovered = [np.sum(cluster == 0)/fullen]
            
        
        ci = np.copy(c)
        nextmeas = []
        numrecover = []
        
        while True:
            cnum = np.sum(latdist<=confcut,axis = 1)
            cnum[~np.isin(latentnames,confiset)] = 0
            argmax = np.argmax(cnum)
            measp = latentnames[argmax]
            addprot = cnum[argmax]
            cprots = latdist[argmax]<=confcut
            
            cluster[np.isin(allnames,latentnames[cprots])] = ci
            nextmeas.append(measp)
            numrecover.append(addprot)
            recovered.append(addprot/fullen)
            
            print ci, recovered[-1], addprot
            
            ci += 1
            latdist = latdist[~cprots][:,~cprots]
            latentnames = latentnames[~cprots]
            confiset = confiset[np.isin(confiset, latentnames)]
            

            if ci == maxnum or len(latentnames) == 0: #recovered[-1] < 0.00025 or numrecover[-1] < 3: , 1259
                break

        recoveredlist.append(recovered)
        np.savetxt('Proposed_next_measurements_cut'+str(confcut)+add+'.txt', np.array([nextmeas,numrecover]).T, fmt = '%s')
        np.savetxt('Recovered_proposed_next_measurement_cut'+str(confcut)+add+'.txt', np.array([allnames,cluster]).T, fmt = '%s')

else:
    recoveredlist = []
    files = sys.argv[1].split(',')
    confcutset = sys.argv[2].split(',')
    add = ''
    for cfile in files:
        clusterfile = np.genfromtxt(cfile,dtype = str)
        clusters, recovered = np.unique(clusterfile[:,-1].astype(int), return_counts = True)
        recovered = recovered[clusters >= 0]/float(len(clusterfile))
        recoveredlist.append(recovered)
        print cfile, 'additional clusters\n', len(recovered), 'for', np.sum(recovered)
        

    
fig = plt.figure(figsize = (4,4))
ax = fig.add_subplot(111)
xticks = []
for r, recovered in enumerate(recoveredlist):
    xva = np.arange(len(recovered), dtype = float)
    ax.plot(xva, np.cumsum(recovered), c = 'k')
    ax.fill_between(xva, 0, np.cumsum(recovered), color = 'grey', alpha = 0.4)
    ax.text(len(recovered), np.sum(recovered)*1.01, str(confcutset[r]), ha = 'right', va = 'bottom')
    xticks.append(len(recovered))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Additionally measured RBPs')
ax.set_ylabel('Percentage Eukaryotic RBPs')
ax.set_yticks(np.linspace(0,1,11))
ax.set_yticklabels([0,10,20,30,40,50,60,70,80,90,100])
ax.grid()
if '--logx' in sys.argv:
    ax.set_xscale('symlog')
    xticks = np.append([0, 1,10,100,1000], xticks)
    xticklabels = xticks.astype(str)
    xticklabels[xticklabels == '0'] = 'Measured'
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, rotation = 60)
ax.set_ylim([0,1])

if '--savefig' in sys.argv:
    outname = 'Additional_measurements_cut'+add+'-'+'-'.join(np.array(confcutset).astype(str))+'.jpg'
    outname = outname.replace(' ', '-')
    print outname
    fig.savefig(outname, bbox_inches = 'tight', dpi = 300)

plt.show()

