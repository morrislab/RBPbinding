import numpy as np
import sys, os
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import glob

# split KH, RRM

def numkmer(x, a, m):
    y = (1./a)*np.log(m*a*x+1.)
    return y

def seqcluster(latdist, confident, domaintype, confcut):
    c = 0
    sclusters = np.zeros(len(latdist))
    for do, dom in enumerate(np.unique(domaintype)):
        pident = np.where(domaintype==dom)[0]
        dlat = latdist[domaintype==dom][:, domaintype==dom]
        dlat = dlat<=confcut
        dconf = confident[domaintype == dom]
        while True:
            neigh = np.sum(dlat, axis = 1)
            if dconf.any():
                neigh[~dconf] = 0
            nclust = dlat[np.argmax(neigh)]
            sclusters[pident[nclust]] =c
            
            dlat = dlat[~nclust][:,~nclust]
            pident = pident[~nclust]
            dconf = dconf[~nclust]
                
            c +=1
            if len(pident) == 0:
                break
    
    nmeas = len(np.unique(sclusters[confident]))
    nc = len(np.unique(sclusters))
            
    return sclusters, nmeas, nc


if __name__ == '__main__':
    if '--analyze' in sys.argv:
        species = np.genfromtxt(sys.argv[1], dtype = str)
        kings = species[:, int(sys.argv[2])]
        species = species[:,0]

        evodistance = np.genfromtxt(sys.argv[3])
        evodistance = np.around(evodistance, 3)
        evospecies = open(sys.argv[3], 'r').readline().strip().replace("'",'').split()[1:]
        sort = []
        for s, spec in enumerate(species):
            sort.append(list(evospecies).index(spec))
        evospecies = np.array(evospecies)[sort]
        evodistance = evodistance[sort][:,sort]

        latentrep = sys.argv[4]
        simcut = float(sys.argv[5])
        disimcut = float(sys.argv[6])


        latentclose = np.genfromtxt(sys.argv[7], dtype = str)
        confcut = float(sys.argv[8])
        latentclose = latentclose[np.argsort(latentclose[:, 0])]
        confident = latentclose[:,-1].astype(float)<=confcut
        latentclose = latentclose[:,0]
        
        protein_simfile = np.load(sys.argv[9])
        protsim = protein_simfile['identmat']
        protsimnames = protein_simfile['names']


        outname = 'Protseqcomparespec_'+os.path.splitext(sys.argv[1])[0]+'_'+os.path.splitext(sys.argv[4])[0]
            
        protdomains = np.zeros((8,len(protsimnames)), dtype = np.int8)
        protdomains2 = np.zeros((8,len(protsimnames)), dtype = np.int8)
        realprotnames = np.copy(protsimnames)
        domaincounts = np.zeros(len(protsimnames))
        
        # generate arrays that determine which proteins can be compared to each other based on their domain composition
        # for example 1 rrm can be compared to 1 and 2 rrms but not 3
        # remove rbps with no info on domains
        remove = np.zeros(len(protsimnames)) == 0
        for p, psname in enumerate(protsimnames):
            psname = psname.split('__')
            protsimnames[p] = psname[0]
            
            if ('KH' not in psname[-1] and 'RRM' not in psname[-1]) or 'ENSG0' in psname[2]:
                remove[p] = False
            
            if '-' in psname[-1]:
                dtypes, dtnum = np.unique(psname[-1].split('-'), return_counts = True)
            else:
                dtypes, dtnum = [psname[-1]], [1]
        
            if 'RRM' in psname[-1]:
                lok = [min(max(0,x),3) for x in range(-2+min(5,dtnum[list(dtypes).index('RRM')]), 1+min(5,dtnum[list(dtypes).index('RRM')]))]
                protdomains[lok,p] = 1
                protdomains2[-1+min(4,dtnum[list(dtypes).index('RRM')]), p] = 1
            if 'KH' in psname[-1]:
                lok = [4+min(max(0,x),3) for x in range(-2+min(5,dtnum[list(dtypes).index('KH')]), 1+min(5,dtnum[list(dtypes).index('KH')]))]
                protdomains[lok,p] = 1
                protdomains2[3+min(4,dtnum[list(dtypes).index('KH')]), p] = 1
            domaincounts[p] = np.sum(dtnum)
            realprotnames[p] = psname[2]+'__'+psname[3]
        
        protsort = np.argsort(protsimnames)
        protsimnames = protsimnames[protsort]
        realprotnames = realprotnames[protsort]
        protsim = protsim[protsort][:, protsort]
        protdomains = protdomains[:, protsort] >=1
        protdomains2 = protdomains2[:, protsort] >=1
        
        confident = confident[np.isin(latentclose, protsimnames)]
        latentclose = latentclose[np.isin(latentclose, protsimnames)]
        
        protsim = protsim[remove][:, remove]
        protsimnames = protsimnames[remove]
        realprotnames = realprotnames[remove]
        protdomains = protdomains[:, remove]
        protdomains2 = protdomains2[:, remove]
        confident = confident[remove]
        latentclose = latentclose[remove]
        print np.array_equal(latentclose, protsimnames)
        
        
        latent = []
        latentnames = []
        latdomains = []
        for s, spec in enumerate(species):
            slat = np.load(spec+'/'+latentrep)
            snames = slat['names']
            ssort = np.argsort(snames)[np.isin(np.sort(snames), protsimnames)]
            sprof = slat['profiles']
            sprof = sprof[ssort]
            snames = snames[ssort]
            latent.append(sprof)
            latentnames.append(snames)
        
        protsort = np.isin(np.sort(protsimnames), np.concatenate(latentnames))
        protsimnames = protsimnames[protsort]
        realprotnames = realprotnames[protsort]
        protsim = protsim[protsort][:, protsort]
        protdomains = protdomains[:, protsort]
        protdomains2 = protdomains2[:, protsort]
        confident = confident[protsort]
        latentclose[protsort]
        
        # remove proteins with short sequences, may be false detections from pfam scanning
        if '--correct_length' in sys.argv:
            lcorfile = sys.argv[sys.argv.index('--correct_length')+1]
            lcorlength = int(sys.argv[sys.argv.index('--correct_length')+2])
            lcorfile = open(lcorfile, 'r').readlines()
            lcrname = []
            lclen = []
            for l, line in enumerate(lcorfile):
                if line[0] == '>':
                    lcrname.append(line[1:].split('__')[0])
                else:
                    lclen.append(len(line.strip()))
            lclen = np.array(lclen)
            lcrname = np.array(lcrname)
            lcrname = lcrname[lclen >= lcorlength]
            
            
            print len(protsimnames)
            for s, spec in enumerate(species):
                sort = np.isin(latentnames[s], lcrname)
                latentnames[s] = latentnames[s][sort]
                latent[s] = latent[s][sort]
            
            sort = np.isin(protsimnames, lcrname)
            protsimnames = protsimnames[sort]
            realprotnames = realprotnames[sort]
            protsim = protsim[sort][:, sort]
            protdomains = protdomains[:, sort]
            protdomains2 = protdomains2[:, sort]
            confident = confident[sort]
            latentclose = latentclose[sort]
            print len(protsimnames)
        
        # generate 3 matrices that mask only rbps that can be compared depending on their domain composition
        protdommat = np.dot(protdomains2.T, protdomains) > 0
        protsim[~protdommat] = 0.
        protdomains = np.array([np.sum(protdomains[:4],axis = 0)>0, np.sum(protdomains[4:],axis = 0)>0])
        
        
        overlap = np.zeros((2, len(species), len(species), 4))
        overlap11 = np.zeros((2, len(species), len(species), 4))
        
        for s, spec in enumerate(species):
            specnames = np.isin(protsimnames, latentnames[s])
            print np.array_equal(latentnames[s], protsimnames[specnames])
            for t, tpec in enumerate(species):
                tpecnames = np.isin(protsimnames, latentnames[t])
                for p, pdom in enumerate(protdomains):
                    # get idmat
                    dommat = protdomains[p]
                    idmat = protsim[specnames*dommat][:, tpecnames*dommat]
                    # get latentdists
                    sprot = np.isin(latentnames[s], protsimnames[specnames*dommat])
                    tprot = np.isin(latentnames[t], protsimnames[tpecnames*dommat])
                    if spec == tpec:
                        overlap11[p,s,t,3] = np.sum(sprot)
                    else:
                        if np.sum(sprot) > 0:
                            if np.sum(tprot) > 0:
                                print spec, tpec
                                latdist = cdist(latent[s][sprot], latent[t][tprot], 'cosine')
                                #print np.shape(latdist)
                                # identify closest rbp in t to protein in s
                                bestst = np.argmax(idmat, axis = 1)
                                # identify closest rbp in s to t
                                bestts = np.argmax(idmat, axis = 0)
                                # identify s proteins with one-to-one orthologs
                                onetoone = bestts[bestst] == np.arange(len(bestst))
                                #onetomanyt = ~np.isin(np.arange(len(bestts)), bestst[onetoone])
                                onetomany = ~onetoone
                                # get identity
                                bestid = idmat[np.arange(len(bestst)), bestst]
                                
                                # check their confidence to measured proteins
                                allconfident = confident[np.isin(latentclose, protsimnames[specnames*dommat])]
                                allconfident = allconfident *(bestid>0)
                                #print np.sum(allconfident)
                                # check binding similarity of closest in s on t
                                similarity = latdist[np.arange(len(latdist), dtype = int), bestst]
                                # count how many are non confident
                                overlap11[p,s,t,0] = np.sum(~allconfident[onetoone]) 
                                overlap[p,s,t,0] = np.sum(~allconfident[onetomany])
                                #print np.sum(~allconfident[onetoone]), np.sum(~allconfident[~onetoone])
                                # count how many are similar
                                overlap11[p,s,t,3] = np.sum(similarity[onetoone*allconfident]<=simcut) 
                                overlap[p,s,t,3] = np.sum(similarity[onetomany*allconfident]<=simcut)
                                #print similarity[onetoone*allconfident], np.sum(similarity[onetoone*allconfident]<=simcut),similarity[onetomany*allconfident], np.sum(similarity[onetomany*allconfident]<=simcut)
                                # count how many are ambiguous
                                overlap11[p,s,t,2] = np.sum((similarity[onetoone*allconfident]>simcut) & (similarity[onetoone*allconfident]<=disimcut) )
                                overlap[p,s,t,2] = np.sum((similarity[onetomany*allconfident]>simcut) & (similarity[onetomany*allconfident]<=disimcut) )
                                # count how many are dissimilar
                                overlap11[p,s,t,1] = np.sum(similarity[onetoone*allconfident]>disimcut)
                                overlap[p,s,t,1] = np.sum(similarity[onetomany*allconfident]>disimcut)
                            else:
                                overlap[p,s,t,0] = np.sum(sprot)
                            
                        #print spec, tpec, overlap[:,s,t,:]
                        #print overlap11[:,s,t,:]
        
        np.savez_compressed(outname+'.npz', overlap11 = overlap11, overlap = overlap, species = species, kings = kings)
        sys.exit()
    elif '--load' in sys.argv:
        npload = np.load(sys.argv[1])
        overlap11 = npload['overlap11']
        overlap = npload['overlap'] 
        species = npload['species']
        if '--outname' in sys.argv:
            outname = sys.argv[sys.argv.index('--outname')+1]
        else:
            outname = os.path.splitext(sys.argv[1])[0]
        speciest = np.genfromtxt(sys.argv[2], dtype = str)
        speciest = speciest[:,0]
        sort = []
        for spt in speciest:
            sort.append(list(species).index(spt))
        
        overlap11sort= overlap11[:,sort]
        overlap11sort= overlap11sort[:,:,sort]
        overlapsort= overlap[:,sort]
        overlapsort= overlapsort[:,:,sort]
        #kings = species[:, int(sys.argv[3])]
        #species = species[:,0]

        evodistance = np.genfromtxt(sys.argv[3])
        evodistance = np.around(evodistance, 3)
        evospecies = open(sys.argv[3], 'r').readline().strip().replace("'",'').split()[1:]
        sort = []
        for s, spec in enumerate(species):
            sort.append(list(evospecies).index(spec))
        evospecies = np.array(evospecies)[sort]
        evodistance = evodistance[sort][:,sort]

    # visualize overlap for one-to-one and other homologs
        # identify similar, ambiguous, and dissimilar and non-identfyable
        # non ident is when too far from training
    # do for kh and rrm
    # plot time resolution for all four
    # average within clade and between clades
    '''
    intspecies = ['Homo_sapiens', 'Schizosaccharomyces_pombe', 'Caenorhabditis_elegans', 'Arabidopsis_thaliana', 'Plasmodium_falciparum'] 
    
    for d, domain in enumerate(['RRM', 'KH']):
        fig2 = plt.figure(figsize = (6.5, 9), dpi = 200)
        ax2a = fig2.add_subplot(221)
        ax2a.set_xscale('symlog', base = 2, linthresh = 100)
        ax2b = fig2.add_subplot(222)
        ax2b.set_xscale('symlog', base = 2, linthresh = 100)
        ax2a2 = fig2.add_subplot(223)
        ax2a2.set_xscale('symlog', base = 2, linthresh = 100)
        ax2b2 = fig2.add_subplot(224)
        ax2b2.set_xscale('symlog', base = 2, linthresh = 100)
        ax2a.spines['top'].set_visible(False)
        ax2a.spines['right'].set_visible(False)
        ax2b.spines['top'].set_visible(False)
        ax2b.spines['right'].set_visible(False)
        ax2a2.spines['top'].set_visible(False)
        ax2a2.spines['right'].set_visible(False)
        ax2b2.spines['top'].set_visible(False)
        ax2b2.spines['right'].set_visible(False)
        for s, speci in enumerate(species):
            if speci in intspecies:
                ax2a.scatter(evodistance[s]+1, overlap11[d][s][:, 3]/overlap11[d][s][s, 3], label = speci)
                ax2b.scatter(evodistance[s]+1, (overlap11[d][s][:, 3]+overlap[d][s][:, 3])/overlap11[d][s][s, 3], label = speci)
                ax2a2.scatter(evodistance[s]+1, overlap11[d][s][:, 1]/overlap11[d][s][s, 3], label = speci)
                ax2b2.scatter(evodistance[s]+1, (overlap11[d][s][:, 1]+overlap[d][s][:, 1])/overlap11[d][s][s, 3], label = speci)
        ax2a.legend(prop={'size' : 6})
        ax2b.legend(prop={'size' : 6})
        ax2a2.legend(prop={'size' : 6})
        ax2b2.legend(prop={'size' : 6})
    '''        
    
    for d, domain in enumerate(['RRM', 'KH']):
        fig = plt.figure(figsize = (len(speciest)*0.25, len(speciest)*0.5), dpi = 50)
        for s, speci in enumerate(speciest):
            ax = fig.add_subplot(len(speciest), 1,s+1)
            ax.bar(np.arange(1,len(speciest)+1), np.sum(overlapsort[d][s], axis = 1), width = 0.4, color = 'silver')
            ax.bar(np.arange(1,len(speciest)+1), np.sum(overlapsort[d][s][:, 1:], axis = 1), width = 0.4, color = 'k')
            ax.bar(np.arange(1,len(speciest)+1), np.sum(overlapsort[d][s][:, 2:], axis = 1), width = 0.4, color = 'grey')
            ax.bar(np.arange(1,len(speciest)+1), np.sum(overlapsort[d][s][:, 3:], axis = 1), width = 0.4, color = 'firebrick')
            ax.bar(np.arange(1,len(speciest)+1)+0.415, np.sum(overlap11sort[d][s], axis = 1), width = 0.4, color = 'silver')
            ax.bar(np.arange(1,len(speciest)+1)+0.415, np.sum(overlap11sort[d][s][:, 1:], axis = 1), width = 0.4, color = 'k')
            ax.bar(np.arange(1,len(speciest)+1)+0.415, np.sum(overlap11sort[d][s][:, 2:], axis = 1), width = 0.4, color = 'grey')
            ax.bar(np.arange(1,len(speciest)+1)+0.415, np.sum(overlap11sort[d][s][:, 3:], axis = 1), width = 0.4, color = 'firebrick')
            ax.text(s+1.615, overlap11sort[d][s][s,3], str(int(overlap11sort[d][s][s,3])), ha = 'left', va = 'top', fontsize = 8)
            #ax.bar([s+1], [overlap[s][::-1][s]], color = 'maroon')
            #ax.bar(np.arange(1,len(species)+1), -underlap[s][::-1], color = 'silver')
            
            ax.plot([0.5,len(species)+0.5], [0,0], lw = 2., color = 'k')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(labelbottom = False, bottom = False)
            ax.set_yticks([0.5*overlap11sort[d][s][s,3]])
            ax.set_yticklabels([speci], rotation = 0)
            ax.set_xlim([0.5,len(speciest)+1.])
            if s == len(speciest)-1:
                ax.tick_params(labelbottom = True)
                ax.set_xticks(np.arange(1,len(speciest)+1))
                ax.set_xticklabels(speciest, rotation = 90)
                
        fig.subplots_adjust(hspace = 0.01)
        if '--savefig' in sys.argv:
            print outname+ domain + '.jpg'
            fig.savefig(outname+ domain + '.jpg', bbox_inches = 'tight', dpi = 250)

        plt.show()

