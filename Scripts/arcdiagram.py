# arcdiagram.py
import numpy as np
import sys, os
import matplotlib.pyplot as plt 
from matplotlib import cm
import pandas as pd
import logomaker as lm

def ic(pw):
    icout = pw*np.log2(pw/0.25)
    icout[icout < -2] = -2
    return icout

def make_arc(a,b):
    x = np.absolute(b-a)/2.*np.sin(np.linspace(-np.pi/2., np.pi,360))+(b+a)/2.
    y = np.sqrt(((b-a)/2.)**2-((2.*x-b-a)/2.)**2)
    return x, y

def arcdiagram(positions, connections, ax, sizes = None, colornodes = None, nodetext = None, colorconnections = None, height = 'equal'):
    if height == 'equal' or height == 'radius':
        ax.set_ylim([-.1, np.amax(np.diff(connections))/1.8])
        if height == 'equal':
            def scale(r):
                return np.amax(np.diff(connections))/2./r
        else:
            def scale(r):
                return 1.
    elif isinstance(height, float):
        ax.set_ylim([-0.1, np.amax(np.diff(connections))/1.8])
        def scale(r):
            return 1./height
    
    ax.set_xlim([-0.5,len(positions)-0.5])
    ax.set_axis_off()
    if sizes is None:
        sizes = np.ones(len(positions))*200
    else:
        sizes = np.sqrt(sizes) * 150
        sizes[sizes > 1050] = 1050
    if colornodes is None:
        colornodes = 'lightgrey'
    if colorconnections is None:
        colorconnections = ['k' for i in range(len(connections))]
    
    
    for c, con in enumerate(connections):
        x, y = make_arc(con[0], con[1])
        y = scale(np.diff(con))*y
        ax.plot(x, y, color = colorconnections[c], lw = 2., zorder = -1)
    ax.scatter(positions, np.zeros(len(positions)), color = colornodes, s = sizes, edgecolor = 'k', zorder = 1)
    if nodetext is not None:
        for n, nodet in enumerate(nodetext):
            if nodet != 0:
                ax.text(n, -0.1, str(nodet), color = 'k', ha = 'center', va = 'bottom', fontsize = min(5.*np.sqrt(nodet), 12))
    
    return 

def plotpwm(pwm, ax, facecolor = None, ylim = [0,2], showy = False, showx = False, frame = True):
    #print pwm
    if pwm is not None:
        pwm = pd.DataFrame({'A':pwm[:,0],'C':pwm[:,1], 'G':pwm[:,2], 'U':pwm[:,3]})        
        pwm = ic(pwm)
        lm.Logo(pwm, ax = ax)
    ax.set_ylim(ylim)
    if frame == False:
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
    if showx:
        ax.spines['bottom'].set_visible(True)
    if showy:
        ax.spines['left'].set_visible(True)
        ax.set_yticks([ylim[0], np.sum(ylim)/2., ylim[1]])
    
    ax.tick_params('both', left=showy, bottom = showx, labelleft = showy, labelbottom = showx)
    if facecolor is not None:
        ax.set_facecolor(facecolor)
    return ax

def readinpwm(pwfile):
    obj = open(pwfile,'r').readlines()
    pws = []
    pwname = []
    for l, line in enumerate(obj):
        line = line.strip().split()
        if len(line) != 0:
            if line[0] == 'Motif':
                if len(pwname) > 0:
                    pws.append(np.array(pw, dtype = float))
                pwname.append(line[1])
                pw = []
            if line[0].isdigit():
                pw.append(line[1:])
    pws.append(np.array(pw, dtype = float))
    return pwname, pws

##### generate common motif for cluster!!!
def clustermot(pwms, mclustercomp):
    clpwms = []
    clmots = []
    for c, clust in enumerate(mclustercomp):
        cpwm = combinepwms(pwms, clust)
        cmot = abbreviatepwm(cpwm)
        cmot = cleanedgesmotif(cmot)
        clpwms.append(cpwm)
        clmots.append(cmot)
    return clpwms, clmots



def cleanedgesmotif(nmotif):
    #print nmotif
    cleanl = 0
    cleanr = 0
    stend = [0,len(nmotif)-1]
    for nn in range(len(nmotif)):
        if cleanl == 0:
            if nmotif[nn] != 'N':
                cleanl += 1
        if cleanl == 1 :
            stend[0] = nn
            cleanl += 1
        if cleanr == 0:
            if nmotif[-nn-1] != 'N':
                cleanr += 1
        if cleanr == 1 :
            stend[1] = len(nmotif)-nn
            cleanr += 1
    #print stend
    nmotif = nmotif[stend[0]:stend[1]]
    #print nmotif
    return nmotif

def kldiv(a1, a2):
    out = np.sum(a1*np.log(a1/a2)+a2*np.log(a2/a1))
    return out

# performs agglomerative clustering of pwms
def combinepwms(pset, clust, weighted = True):
    psetclust = []
    clust, psetweight = np.unique(clust, return_counts = True)
    for clu in clust:
        psetclust.append(pset[clu])
    if weighted == False:
        psetclust = np.ones(len(clust))
    while len(psetclust) > 1:
        klds = []
        klpwms = []
        pairs = []
        for c, clu in enumerate(psetclust):
            for d in range(c+1, len(psetclust)):
                clw = psetclust[d]
                kld, klofs = alignpwms(clu, clw)
                pairs.append([c,d])
                klds.append(kld)
                #### klofs were wrong in alignpwms
                klpwms.append(klofs)
        bekl = np.argmin(klds)
        npsclust = []
        for c, clu in enumerate(psetclust):
            if c not in pairs[bekl]:
                npsclust.append(clu)

        order = np.argsort(psetweight[pairs[bekl]])
        offset = int(klpwms[bekl]*np.diff(order))
        npsclust.append(meanpwm(psetclust[pairs[bekl][order[0]]], psetclust[pairs[bekl][order[1]]], offset, psetweight[pairs[bekl]][order]))
        #print npsclust[-1]
        #print meanpwm(psetclust[pairs[bekl][0]], psetclust[pairs[bekl][1]], klpwms[bekl], psetweight[pairs[bekl]])
        if weighted:
            #print psetweight[~np.isin(np.arange(len(psetweight)), pairs[bekl])], np.sum(psetweight[pairs[bekl]])
            psetweight = np.append(psetweight[~np.isin(np.arange(len(psetweight)), pairs[bekl])], [np.sum(psetweight[pairs[bekl]])])
        else:
            psetweight = np.ones(len(npsclust))
        psetclust = npsclust
        
    return psetclust[0]
            
    
def meanpwm(pa, pb, ofs, weights):
    la = len(pa)
    lb = len(pb)
    pfin = np.zeros((len(pb),4))
    pfin += pb* weights[1]
    pfin[min(0,-(lb+ofs)):min(lb,la-ofs)] += pa[max(0,ofs):min(la,ofs+lb)]*weights[0]
    pnorm = np.ones(lb)*weights[1]
    pnorm[min(0,-(lb+ofs)):min(lb,la-ofs)] +=weights[0]
    return (pfin.T/pnorm).T

def alignpwms(x, y):
    lenx = len(x)
    leny = len(y)
    if lenx > leny:
        a = x
        b = y
        offdir = 1
    else:
        a = y
        b = x
        #### We need to change the offset if we switch the order.
        offdir = -1
    #print a
    #print b
    klset = []
    offset = []
    for i in range(len(a)-len(b)+1):
        pcheck1 = np.copy(b)
        pcheck2 = np.ones((len(b), 4))*0.25
        pcheck2 = a[i:i+len(b)]
        klset.append(kldiv(pcheck1, pcheck2))
        offset.append(i)
    for i in range(2,len(b)):
        pcheck1 = np.copy(b)
        pcheck2 = np.ones((len(b), 4))*0.25
        pcheck2[-i:] = a[:i]
        klset.append(kldiv(pcheck1, pcheck2))
        offset.append(i-len(b))
        pcheck3 = np.ones((len(b), 4))*0.25
        pcheck3[:i] = a[-i:]
        klset.append(kldiv(pcheck1, pcheck3))
        offset.append(len(a)-i)
    bl = np.argmin(klset)
    bestoff = offset[bl]*offdir
    bestkl = klset[bl]
    return bestkl, bestoff
        
        

def abbreviatepwm(pwm):
    motabb = ''
    for pwnum, pwpos in enumerate(pwm):
        sortpw = np.argsort(-pwpos)
        abbvec = np.zeros(4)
        if pwpos[sortpw[0]] >= 0.47:
            if pwpos[sortpw[1]] < 0.8*0.47:
                abbvec[sortpw[0]] = 1
            else:
                abbvec[sortpw[:2]] = 1
        elif pwpos[sortpw[0]] >= 0.3:
            if pwpos[sortpw[1]] < 0.8*.3:
                abbvec[:] = 0
            elif pwpos[sortpw[1]] >= 0.8*0.3:
                if pwpos[sortpw[2]] < 0.8*0.3:
                    abbvec[sortpw[:2]] = 1
                elif pwpos[sortpw[2]] >= 0.8*0.3:
                    abbvec[sortpw[:3]] = 1
        abbrev = (nucabbmat == abbvec).all(axis = 1).nonzero()[0]
        #print abbrev, abbvec
        motabb += nucleotideabbreviations[abbrev[0]]
    return motabb



    
'''
positions = np.arange(10)
connections = [[1,2],[1,4], [2,9]]
sizes = np.random.random(10)*30.
colors = cm.rainbow(np.linspace(0,1,10))
colorconnections = cm.rainbow(np.random.random(3))
'''

if __name__ == '__main__':
    species = open(sys.argv[1], 'r').readline().strip().split()[14:]
    
    clusterprot = np.genfromtxt(sys.argv[1], dtype = str)
    clusterprot = clusterprot[:,1:]
    clusterrealname = clusterprot[:,0]
    clusterid = clusterprot[:, 1].astype(int)
    clusterassign, clusterindex = np.unique(clusterid, return_index = True)
    clustermeasbool = clusterprot[clusterindex,2] == '1'
    
    outname = os.path.splitext(sys.argv[1])[0]
    cluster = clusterprot[:,12:].astype(int)[clusterindex]
    #clustercenter = clusterprot[:,3]
    closestdist = clusterprot[:,3].astype(float)
    closestmeas = clusterprot[:,4]
    clustermeasured = cluster[clustermeasbool]
    
    clusterlen = clusterprot[:,5].astype(int)
    clusterage = clusterprot[:,6].astype(float)[clusterindex]
    seqage = clusterprot[:,7].astype(float)[clusterindex]
    family = clusterprot[:,8].astype(int)[clusterindex]
    evolutionids = clusterprot[:,9].astype(float)
    evolutiontype = clusterprot[:,10].astype(int)
    pspecies = clusterprot[:,11].astype(int)
    evolutiontail = evolutiontype[clusterindex]
    evolutionid = evolutionids[clusterindex]
    

    
    pwmfile = sys.argv[2]
    pwmnames, pwms = readinpwm(pwmfile)

    def alpharatio(string):
        cla = 0.
        for l in string:
            if l.isalpha():
                cla +=1.
        return cla/float(len(string))
    
    speciesrank = ['Homo_sapiens', 'Mus_musculus', 'Macaca_mulatta', 'Gallus_gallus', 'Danio_rerio', 'Xenopus_tropicalis', 'Saccharomyces_cerevisiae', 'Schizosaccharomyces_pombe', 'Drosophila_melanogaster', 'Drosophila_ananassae', 'Caenorhabditis_elegans', 'Caenorhabditis_briggsae', 'Arabidopsis_thaliana', 'Arabidopsis_lyrata', 'Zea_mays', 'Tetraodon_nigroviridis', 'Oryza_sativa', 'Canis_familiaris', 'Latimeria_chalumnae', 'Geospiza_fortis', 'Musca_domestica']
    speciesdic = {}
    for s, spec in enumerate(species):
        if spec in speciesrank:
            speciesdic[s] = speciesrank.index(spec)
        else:
            speciesdic[s] = len(species)
    
    # generate mixture of pwms for cluster from each closest measured protein 
    clusterpwms = []
    clusterpresent = []
    clustername = []
    all_names = []
    stringlen = np.vectorize(len)
    alphacontent = np.vectorize(alpharatio)
    for c, cla in enumerate(clusterassign):
        clnames = [cln.rsplit('__',1)[0] for cln in clusterrealname[clusterid == cla]]
        cldomains = [cln.rsplit('__',1)[1] for cln in clusterrealname[clusterid == cla]]
        #clnames, clncount = np.unique(np.char.upper(clnames), return_counts = True)
        all_names.append(clnames)
        cldomains, cldnum = np.unique(cldomains, return_counts = True)
        ranking = [speciesdic[ps] for ps in pspecies[clusterid == cla]]
        clustername.append(clnames[np.lexsort((stringlen(clnames), ranking))[0]]+'['+cldomains[np.argmax(cldnum)]+']')
        clpwms = []
        clpresent = []
        clustermask = clusterid == cla
        if '--meanPWM' in sys.argv or '--showall' in sys.argv:
            allmeasured = [pwmnames.index(clc) for clc in closestmeas[clustermask]]
            pwmmean = combinepwms(pwms, allmeasured)
        for s, spec in enumerate(species):
            clusterspecies = pspecies[clustermask] == s
            if np.sum(clusterspecies)>0:
                if '--meanPWM' in sys.argv:
                    clpwms.append(pwmmean)
                else:
                    measprot = closestmeas[clustermask][clusterspecies]
                    measprot, meanum = np.unique(measprot, return_counts = True)
                    clpwms.append(pwms[pwmnames.index(measprot[np.argmax(meanum)])])
                clcoin = closestdist[clustermask][clusterspecies] < 0.2
                clpresent.append(clcoin.any())
            else:
                clpwms.append(None)
                clpresent.append(False)
        clusterpwms.append(clpwms)
        clusterpresent.append(clpresent)
        '''
        if we want to assign the cluster a mixture, small differences might get lost
        if clustermeasbool[c]:
            allmeasured = [pwmnames.index(clc) for clc in clustercenter[clusterid == cla]]
            clusterpwms.append(combinepwms(pwms, allmeasured))
        else:
            clusterpwms.append(None)
        '''
    clustername = np.array(clustername)
    
    # assign cluster families
    clusterfamily = [[]]
    clfam = []
    assign = True
    while True:
        lenmax = len(clfam)
        if assign:
            if len(clusterfamily) > 0 and lenmax > 0:
                clusterfamily.append(clfam)
            lenmax = 0
            assign = False
            notinfam = np.setdiff1d(clusterassign, np.concatenate(clusterfamily))
            if len(notinfam) == 0:
                break
            clfam = notinfam[[0]]
        
        # check for every cluster that has clfam as parent and check for every cluster that is parent to any of clfam clusters
        connected = np.append(family[clfam], clusterassign[np.isin(family,clfam)])
        clfam = np.unique(np.append(clfam, connected))
        clfam = clfam[clfam!=-1]
        #print lenmax, len(clfam), clfam
        if len(clfam) == lenmax:
            assign = True
    # remove empty list
    clusterfamily = clusterfamily[1:]
    uniq, uninum = np.unique(np.concatenate(clusterfamily), return_counts= True)
    for cf in clusterfamily:
        print len(cf), np.unique(clusterage[cf], return_counts = True), np.unique(seqage[cf], return_counts = True)
        sort = np.lexsort([seqage[cf], clusterage[cf]])
        for s in sort:
            print clustername[cf][s], clusterage[cf][s], '<--', seqage[cf][s], clustername[family[cf][s]], clusterage[family[cf][s]]
    print uninum[uninum>1], uniq[uninum>1]
    print len(np.concatenate(clusterfamily)), len(np.unique(np.concatenate(clusterfamily)))
    # for each category, give out the unique parents the number of times they serve, and the clusterages they serve for
    evolutiontail, family
    for e, evt in enumerate(np.unique(evolutiontail)):
        evtmask = np.where(evolutiontail == evt)[0]
        evtfam = np.unique(family[evolutiontail == evt])
        print '\nCASE', evt
        for ef in evtfam:
            mask = evtmask[family[evtmask] == ef]
            if len(mask) > 5:
                print '\nFam', clustername[ef], clusterage[ef], len(mask), ','.join(cluster[ef].astype(str)), np.sum(cluster[ef])
                for m in mask:
                    print clustername[m], clusterage[m], seqage[m]
    
    #sys.exit()
            
    if '--showall' in sys.argv:
        outname += '_superfamily'
    if '--measuredonly' in sys.argv:
        outname += '_measuredonly'
                

    for c, clfam in enumerate(clusterfamily):
        print c, len(clfam)
        if '--showall' not in sys.argv or '--measuredonly' in sys.argv:
            # remove clusters that don't server as parent if they are not measured
            clfam = np.union1d(clfam[clustermeasbool[clfam]] , family[clfam[clustermeasbool[clfam]]])
            clfam = clfam[clfam>=0]
        print len(clfam)
        if len(clfam)>1:
            
            positions = np.arange(len(clfam))
            orderfam = np.lexsort((seqage[clfam], clusterage[clfam]))[::-1]
            clfam = clfam[orderfam]
            
            #connections = zip(clfam,family[clfam])
            concolors = cm.tab20([12, 10, 6, 8, 15])
            connections = []
            colorconnections = []
            for f, fm in enumerate(family[clfam]):
                if fm in clfam:
                    connections.append([list(clfam).index(fm), f])
                    colorconnections.append(concolors[evolutiontail[clfam[f]]])
            connections = np.array(connections)
            famc, famcount = np.unique(connections[:,0], return_counts = True)
            outcount = [famcount[list(famc).index(x)] if x in famc else 0 for x in positions]
            connections = np.sort(connections, axis = 1)
            rmax = np.amax(np.diff(connections))/2.
            fig = plt.figure(figsize=(len(clfam)*0.5, (len(species)+rmax)*0.25), dpi = 100)
            
            axarc = fig.add_subplot(121)
            yhigh = 0.8 * rmax/(float(len(species))+rmax)
            axarc.set_position([0.1,0.9-yhigh,0.8,yhigh])
            arcdiagram(positions, connections, axarc, sizes = outcount, colornodes = None, nodetext = outcount, colorconnections = colorconnections, height = 'radius')
            
            if '--showall' in sys.argv or '--measuredonly' in sys.argv:
                for d, cl in enumerate(clfam):
                    ax = fig.add_subplot(len(species)+1, len(clfam), d+1)
                    pheight = (0.8-yhigh)/float(len(species)+1)
                    pwidth = 0.8/float(len(clfam))
                    ax.set_position([0.1+d*pwidth,0.9-yhigh-pheight, pwidth, pheight])
                
                    pwm = clusterpwms[cl][np.where(cluster[cl]>0)[0][0]]
                    if clustermeasbool[cl]:
                        plotpwm(pwm, ax, facecolor = 'white')
                    else:
                        plotpwm(pwm, ax, facecolor = 'grey')
                    
                ax2 = fig.add_subplot(1, len(clfam), d+1)
                ax2.set_position([0.1 ,0.1, len(clfam)*pwidth, len(species)*pheight])
                
                ax2.imshow(cluster[clfam].T, cmap = 'Greys', vmin = 0, vmax = 5, aspect='auto')
                
                xlabels = []
                for d, cl in enumerate(clfam):
                    xlabels.append('('+str(cl)+')'+clustername[cl]+' ('+str(clusterage[cl])+' MyA)')
                    for s, spec in enumerate(species):
                        if cluster[cl,s] > 3:
                            tcol = 'white'
                        else:
                            tcol = 'k'
                        if cluster[cl,s] > 0:
                            ax2.text(d, s, str(int(cluster[cl,s])), color = tcol, va = 'center', ha = 'center')
                ax2.set_xticks(np.arange(len(clfam)))
                ax2.set_xticks(np.arange(len(clfam))-0.5, minor = True)
                ax2.tick_params(labelleft = False,left = False, which = 'minor', labelbottom = False, bottom = False)
                ax2.tick_params(left = False, which = 'major', bottom = False)
                ax2.set_xticklabels(xlabels, rotation = 90, va = 'top', ha = 'center')
                
                ax2.set_yticks(np.arange(len(species)))
                ax2.set_yticks(np.arange(len(species))-0.5, minor = True)
                ax2.set_yticklabels(species)
                ax2.grid(which = 'minor', color = 'k')
            else:
                for d, cl in enumerate(clfam):
                    #print d, cl, clusterpwms[cl], len(clfam), cluster[cl]
                    for s, spec in enumerate(species):
                        ax = fig.add_subplot(len(species), len(clfam), d*len(species)+s+1)
                        pheight = (0.8-yhigh)/float(len(species))
                        pwidth = 0.8/float(len(clfam))
                        #print [0.1+d*pwidth,0.9-yhigh-(s+1.)*pheight,pwidth,pheight]
                        ax.set_position([0.1+d*pwidth,0.9-yhigh-(s+1.)*pheight,pwidth,pheight])
                        if clusterpwms[cl][s] is not None:
                            if clusterpresent[cl][s]:
                                plotpwm(clusterpwms[cl][s], ax)                
                            else:
                                if clustermeasbool[cl]:
                                    plotpwm(clusterpwms[cl][s], ax, facecolor = 'lightgrey')
                                else:
                                    plotpwm(clusterpwms[cl][s], ax, facecolor = 'grey')
                        else:
                            plotpwm(None, ax, facecolor = 'white')
                        
                        if d == 0:
                            ax.tick_params(labelleft = True)
                            ax.set_yticks([1])
                            ax.set_yticklabels([spec], rotation = 0, va = 'center', ha = 'right')
                        if s == len(species)-1:
                            ax.tick_params(labelbottom = True)
                            ax.set_xticks([0])
                            ax.set_xticklabels([''])
                            ax.set_xlabel('('+str(cl)+')'+clustername[cl]+' ('+str(clusterage[cl])+' MyA)', rotation = 90, va = 'top', ha = 'center')
                
                
            if '--savefig' in sys.argv:
                if len(sys.argv) > sys.argv.index('--savefig')+1:
                    if sys.argv[sys.argv.index('--savefig')+1] == 'png':
                        fmt = '.png'
                    else:
                        fmt = '.jpg'
                fig.savefig(outname+'_clusterfamily_arc'+str(c)+fmt, dpi = 150, bbox_inches = 'tight')
                print outname+'_clusterfamily_arc'+str(c)+fmt
            else:
                plt.show()
            plt.close()
    
    
    
    
    
    
    






