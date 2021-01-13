#make_newick_phylo.py
import sys, os
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from functools import reduce
import logomaker as lm
import pandas as pd
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62
from scipy.stats import pearsonr
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.patches import Rectangle



# Read in different files
def readinclustermat(clustermat):
    clustnames = open(clustermat, 'r').readline().strip().split()[1:]
    if '||' in clustnames[0]:
        for c, cname in enumerate(clustnames):
            cname = cname.split('||')
            clustnames[c] = cname[-1]
    clusteringmat = np.genfromtxt(clustermat)
    return np.array(clustnames), clusteringmat

def keepfunc(names, matrix, rightnames):
    keep = []
    for rightname in rightnames:
        keep.append(list(names).index(rightname))
    names = np.array(names)[keep]
    matrix = matrix[keep]
    matrix = matrix[:,keep]
    return names, matrix

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

def readinfasta(ffile):
    obj = open(ffile, 'r').readlines()
    expname = []
    sequence = []
    for l, line in enumerate(obj):
        if line[0] == '>':
            expname.append(line.strip().split('||')[-1])
        else:
            sequence.append(line.strip())
    return expname, sequence


### Plot PWMs
def ic(pw):
    icout = pw*np.log2(pw/0.25)
    icout[icout < -2] = -2
    return icout

def plotpwm(pwm, ax):
    #print pwm
    pwm = pd.DataFrame({'A':pwm[:,0],'C':pwm[:,1], 'G':pwm[:,2], 'U':pwm[:,3]})        
    pwm = ic(pwm)
    lm.Logo(pwm, ax = ax)
    #ax.set_yticks([0,1,2])
    #ax.set_yticklabels([0,1,2])
    ax.set_ylim([0,2])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params('both', left=False, bottom = False, labelleft = False, labelbottom = False)
    return ax








##### generate common PWM for cluster!!!
def clustermot(pwms, idnames, mclustercomp):
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

def combinepwms(pset, clust):
    psetclust = []
    for clu in clust:
        psetclust.append(pset[clu])
    while len(psetclust) > 1:
        klds = []
        klpwms = []
        pairs = []
        for c, clu in enumerate(psetclust):
            for d in range(c+1, len(psetclust)):
                clw = psetclust[d]
                kld, klpwm = alignpwms(clu, clw)
                pairs.append([c,d])
                klds.append(kld)
                klpwms.append(klpwm)
        bekl = np.argmin(klds)
        npsclust = []
        for c, clu in enumerate(psetclust):
            if c not in pairs[bekl]:
                npsclust.append(clu)
        npsclust.append(klpwms[bekl])
        psetclust = npsclust
    return psetclust[0]
            
    
def meanpwm(pa, pb, ofs):
    la = len(pa)
    lb = len(pb)
    pfin = np.zeros((len(pb),4))
    pfin += pb
    pfin[min(0,-(lb+ofs)):min(lb,la-ofs)] += pa[max(0,ofs):min(la,ofs+lb)]
    pnorm = np.ones(lb)
    pnorm[min(0,-(lb+ofs)):min(lb,la-ofs)] +=1.
    return (pfin.T/pnorm).T

def alignpwms(x, y):
    lenx = len(x)
    leny = len(y)
    if lenx > leny:
        a = x
        b = y
    else:
        a = y
        b = x
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
    bestoff = offset[bl]
    bestkl = klset[bl]
    mpwm = meanpwm(a, b, bestoff)
    return bestkl, mpwm
        
nucleotideabbreviations = ['A','C','G','U','M','R','W','S','Y','K','V','H','D','B','N']
nucabbmat = np.array([[1,0,0,0],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,1],
            [1,1,0,0],
            [1,0,1,0],
            [1,0,0,1],
            [0,1,1,0],
            [0,1,0,1],
            [0,0,1,1],
            [1,1,1,0],
            [1,1,0,1],
            [1,0,1,1],
            [0,1,1,1],
            [0,0,0,0]])
    

# function to transform PFM into IUPAC 
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


### assigns protein names to clusters
def clusterclean(measuredname, proteinnames, assigncl):
    assignpn = []
    assignspecies = []
    for m, mn in enumerate(measuredname):
        ap = list(proteinnames[:,0]).index(mn)
        apname = proteinnames[ap,1].replace('_CONSTRUCT', '')+' ('+r"$\it{}$".format(proteinnames[ap,2][0]+'. '+proteinnames[ap,2].split('_')[1])+')'
        assignpn.append(apname)
        assignspecies.append(proteinnames[ap, 2])
    
    assigncl = np.array(assigncl)
    assignpn = np.array(assignpn)
    assignspecies = np.array(assignspecies)
    outname = []
    outpname = []
    for clprt in assigncl:
        if len(clprt) > 1:
            specishort = assignspecies[clprt]
            cprt = []
            addnum = [[]]
            if 'Homo_sapiens' in specishort:
                cprt.append(np.unique(assignpn[clprt[specishort == 'Homo_sapiens']]))
                addnum.append(clprt[specishort == 'Homo_sapiens'])
            if 'Mus_musculus' in specishort:
                cprt.append(np.unique(assignpn[clprt[specishort == 'Mus_musculus']]))
                addnum.append(clprt[specishort == 'Mus_musculus'])
            if 'Drosophila_melanogaster' in specishort:
                cprt.append(np.unique(assignpn[clprt[specishort == 'Drosophila_melanogaster']]))
                addnum.append(clprt[specishort == 'Drosophila_melanogaster'])
            if 'Caenorhabditis_elegans' in specishort:
                cprt.append(np.unique(assignpn[clprt[specishort == 'Caenorhabditis_elegans']]))
                addnum.append(clprt[specishort == 'Caenorhabditis_elegans'])
            if 'Saccharomyces_cerevisiae' in specishort:
                cprt.append(np.unique(assignpn[clprt[specishort == 'Saccharomyces_cerevisiae']]))
                addnum.append(clprt[specishort == 'Saccharomyces_cerevisiae'])
                
            if len(cprt) < 3:
                addnum = np.concatenate(addnum)
                cprt.append(np.unique(assignpn[clprt[~np.isin(clprt, addnum)][:3-len(cprt)]]))
            cprt = np.concatenate(cprt)
            if len(cprt) > 3:
                cprt = np.append(cprt[:3], ['...'])
            
            outname.append(','.join(cprt))

        else:
            outname.append(assignpn[clprt[0]])
    return np.array(outname)


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# aligns protein sequence, calculates identities of different regions (domain, linker) and returns an array that represents domains and identity which can be plotted
def align(s1, s2, l1, l2):
    # check l1 and l2 entries for redundant start and end points
    print l1
    for l in l1:
        if len(np.where(l1 == l)[0])>1:
            l1[np.where(l1 ==l)[0][-1]] += 1
    print l1
    print l2
    for l in l2:
        if len(np.where(l2 == l)[0])>1:
            l2[l2 ==l][-1] += 1
    print l2
        
    alignment = pairwise2.align.localds(s1, s2, matrix, -11, -1)[0]
    lseq1 = len(s1)-1
    lseq2 = len(s2)-1
    alignmat = []
    c1 = 0
    c2 = 0
    al1 = alignment[0]
    al2 = alignment[1]
    indom1 = False
    indom2 = False
    pids1 = []
    pids2 = []
    ls1 = []
    ls2 = []
    pid1 = 0.
    pl1 = 0.
    part1 = 0
    pid2 = 0.
    pl2 = 0.
    part2 = 0
    inseq1 = False
    inseq2 = False
    for a in range(len(al1)):
        #print al1[a], al2[a]
        if al1[a] != '-':
            cha = True
            c1 += 1
            pl1 += 1.
        else:
            cha = False
        if al2[a] != '-':
            pl2 += 1.
            chb = True
            c2 += 1
        else:
            chb = False
        
        if al1[a] == al2[a]:
            pid1 += 1.
            pid2 += 1.
                
        if (c1 == 1 and cha) and (c2 == 1 and chb):
            alignmat.append([-1,-1,-1,-1,-1,-1])
            inseq1 = True
            inseq2 = True
            part1 +=1
            part2 +=1
        elif (c1 == 1 and cha):
            inseq1 = True
            part1 +=1
            if inseq2 == False:
                alignmat.append([-1,-1,-1,0,0,0])
            elif indom2:
                alignmat.append([-1,-1,-1,part2, part2, part2])
            else:
                alignmat.append([-1,-1,-1,part2, 0,0])
        elif (c2 == 1 and chb):
            part2 += 1
            inseq2 = True
            if inseq1 == False:
                alignmat.append([0,0,0,-1,-1,-1])
            elif indom1:
                alignmat.append([part1, part1, part1,-1,-1,-1])
            else:
                alignmat.append([0,0,part1,-1,-1,-1])
        
        if (c1 == lseq1 and cha) and (c2 == lseq2 and chb):
            inseq1 = False
            inseq2 = False
            alignmat.append([-1,-1,-1,-1,-1,-1])
        elif (c1 == lseq1 and cha):
            inseq1 = False
            if inseq2 == False:
                alignmat.append([-1,-1,-1,0,0,0])
            elif indom2:
                alignmat.append([-1,-1,-1,part2, part2, part2])
            else:
                alignmat.append([-1,-1,-1,part2, 0,0])
        elif (c2 == lseq2 and chb):
            inseq2 = False
            if inseq1 == False:
                alignmat.append([0,0,0,-1,-1,-1])
            elif indom1:
                alignmat.append([part1, part1, part1,-1,-1,-1])
            else:
                alignmat.append([0,0,part1,-1,-1,-1])           
        
        elif (c1 in l1 and cha) and (c2 in l2 and chb):
            alignmat.append([-1,-1,-1,-1,-1,-1])
            alignmat.append([-1,-1,-1,-1,-1,-1])
            indom1 = ~indom1
            indom2 = ~indom2
            part1 += 1
            part2 += 1
            pids1.append(pid1/pl1)
            pids2.append(pid2/pl2)
            ls1.append(pl1)
            ls2.append(pl2)
            pid1 = 1.
            pl1 = 1.
            pid2 = 1.
            pl2 = 1.
        elif (c1 in l1 and cha):
            if inseq2 == False:
                alignmat.append([-1,-1,-1,0,0,0])
                alignmat.append([-1,-1,-1,0,0,0])
            elif indom2:
                alignmat.append([-1,-1,-1,part2,part2,part2])
                alignmat.append([-1,-1,-1,part2,part2,part2])
            else:
                alignmat.append([-1,-1,-1,part2,0,0])
                alignmat.append([-1,-1,-1,part2,0,0])
            indom1 = ~indom1
            part1 += 1
            pids1.append(pid1/pl1)
            ls1.append(pl1)
            pid1 = 1.
            pl1 = 1.
        elif (c2 in l2 and chb):
            if inseq1 == False:
               alignmat.append([0,0,0,-1,-1,-1])
               alignmat.append([0,0,0,-1,-1,-1])
            if indom1:
                alignmat.append([part1,part1,part1,-1,-1,-1])
                alignmat.append([part1,part1,part1,-1,-1,-1])
            else:
                alignmat.append([0,0,part1,-1,-1,-1])
                alignmat.append([0,0,part1,-1,-1,-1])
            indom2 = ~indom2
            part2 += 1
            pids2.append(pid2/pl2)
            ls2.append(pl2)
            pid2 = 1.
            pl2 = 1.

        else:
            if inseq1 == False:
                appen1 = [0,0,0]
            elif indom1:
                appen1 = [part1, part1, part1]
            else:
                appen1 = [0, 0, part1]
            
            if inseq2 == False:
                appen2 = [0,0,0]            
            elif indom2:
                appen2 = [part2, part2, part2]
            else:
                appen2 = [part2, 0, 0]
            appen = appen1 + appen2
            alignmat.append(appen)
        #print indom1, indom2, alignmat[-1]
    
    part1 +=1
    pids1.append(pid1/pl1)
    ls1.append(pl1)
    ls2.append(pl2)
    part2 +=1
    pids2.append(pid2/pl2)
    alignmat = np.array(alignmat, dtype = float).T
    alignticks = [[],[]]
    alignticklabels = [[],[]]
    for p in range(1,part1):
        if ls1[p-1] >= 15:
            alignticklabels[0].append(str(int(pids1[p-1]*100))+'%')
            alignticks[0].append(np.median(np.where(alignmat[2] == p)[0]))
        alignmat[:3][alignmat[:3] == p] = 0.1+pids1[p-1]
    for p in range(1,part2):
        if ls2[p-1] >= 15:
            alignticklabels[1].append(str(int(pids2[p-1]*100))+'%')
            alignticks[1].append(np.median(np.where(alignmat[3] == p)[0]))
        alignmat[3:][alignmat[3:] == p] = 0.1+pids2[p-1]
        
    alignmat[alignmat == -1] = 1.1
    
    return alignmat, alignticks, alignticklabels
    

def topz(x,y, t = 100):
    sorty = np.argsort(y)[-t:]
    sortx = np.argsort(x)[-t:]
    return float(len(np.intersect1d(sorty, sortx)))/float(t) #len(np.union1d(sorty, sortx)))

    



    
if __name__ == '__main__':
    
    # provide file for clustering: Motif similarity or sequence identity
    simfile = sys.argv[1]
    idnames, idmat = readinclustermat(simfile)
    idmatmax = np.amax(idmat)
    np.fill_diagonal(idmat, np.amax(idmat))
    idmat = 1. - idmat/np.amax(idmat)
    
    # define clustering parameters
    clusttype = sys.argv[2] # ward, single complete, median, average weighted centroid
    cutoff = sys.argv[3]
    clustername = sys.argv[4]
    clusterticklabels = sys.argv[5].split(',')
    if is_number(cutoff):
        clustercut = float(cutoff)
        clustercut = 1. - clustercut/idmatmax
    else:
        clustercut = None
    clusterticks = 1. - np.array(clusterticklabels, dtype = float)/idmatmax
    
    outname = os.path.splitext(os.path.split(simfile)[1])[0] +'_cluster'+clusttype+'_cutoff'+cutoff
    if '--savefig' in sys.argv:
        outname = sys.argv[sys.argv.index('--savefig')+1]+outname
        dpi = int(sys.argv[sys.argv.index('--savefig')+2])
        
    
    
        
    
    
    # provide protein information
    masterfile = open(sys.argv[6], 'r').readlines()
    masterheader = masterfile[0].strip().split('\t')
    masterinfo = []
    for l, line in enumerate(masterfile):
        if l > 0:
            line = line.strip().split('\t')
            if ',' in line[masterheader.index('Motif_ID')]:
                rncmpt = line[masterheader.index('Motif_ID')].split(',')
            else:
                rncmpt = [line[masterheader.index('Motif_ID')]]
            for r, rnc in enumerate(rncmpt):
                nline = []
                for e, entry in enumerate(line):
                    
                    if ',' in entry:
                        entry = entry.split(',')
                        if entry[0].isdigit():
                            entry = np.array(entry,dtype = int)
                    elif e == 8:
                        entry = [entry]
                    nline.append(entry)
                nline[1] = rnc
                for s in range(7):
                    nline[-s-1] = nline[-s-1][r]
                masterinfo.append(nline)
    masterinfo = np.array(masterinfo)
    
    #provide PWMs sort to idnames
    motiffile = sys.argv[7]
    pwmnames, pwms = readinpwm(motiffile)


    # choose set of RNCMPT experiments
    if '--proteinlist' in sys.argv:
        proteinlist = sys.argv[sys.argv.index('--proteinlist')+1]
        outname += '_plist'+os.path.splitext(os.path.split(proteinlist)[1])[0]
        prots = np.genfromtxt(proteinlist, dtype = str)
        if len(np.shape(prots)) == 2:
            prots = prots[:,1]
        prots = np.sort(prots)
        idnames, idmat = keepfunc(idnames, idmat, prots)
        prots = idnames
    elif '--filterlist' in sys.argv:
        filtercut = 1. - float(sys.argv[sys.argv.index('--filterlist')+1])/idmatmax
        Z = linkage(idmat[np.triu_indices(len(idmat),1)], 'complete')
        mclustercomp, Z = cleanz(Z, filtercut, idnames)
        keep = []
        for m, mclust in enumerate(mclustercomp):
            keep.append(mclust[-1])
        idnames = idnames[keep]
        idmat = idmat[keep][:, keep]
        prots = idnames
    else:
        prots = idnames


    # clustering here
    Z = linkage(idmat[np.triu_indices(len(idmat),1)], clusttype)    
    if clustercut is None:
        mclustercomp = np.arange(len(prots)).reshape(-1,1)
    else:
        mclustercomp, Z = cleanz(Z, clustercut, prots)
    
    mclustlen = []
    for m, mc in enumerate(mclustercomp):
        mclustercomp[m] = np.array(mc,dtype = int)
        mclustlen.append(len(mc))
    mclustlen = np.array(mclustlen)
    
    # remove clusters that are smaller than defines size
    if '--minclustersize' in sys.argv:
        mincsize = int(sys.argv[sys.argv.index('--minclustersize')+1])
        outname +='_minclsize'+str(mincsize)
        cmask = mclustlen >= mincsize
        mclustercomp = list(np.array(mclustercomp)[cmask])
        mclustlen = mclustlen[cmask]
        restprots = np.sort(np.concatenate(mclustercomp)).astype(int)
        idmat = idmat[restprots][:,restprots]
        prots = prots[restprots]
        idnames = idnames[restprots]
        Z = linkage(idmat[np.triu_indices(len(idmat),1)], clusttype)
        mclustercomp, Z = cleanz(Z, clustercut, prots)
        mclustlen = []
        for m, mc in enumerate(mclustercomp):
            mclustercomp[m] = np.array(mc, dtype = int)
            mclustlen.append(len(mc))
        mclustlen = np.array(mclustlen)
    else:
        mincsize = 0


    # choose format of plot
    if '--showsimilarity' in sys.argv:
        similarityfile = sys.argv[sys.argv.index('--showsimilarity')+1]
        outname += '_similaritymat'+os.path.splitext(os.path.split(similarityfile)[1])[0]
        simnames, simmat = readinclustermat(similarityfile)
        nsimmat = np.zeros((len(mclustercomp), len(mclustercomp)))
        for m, mclust in enumerate(mclustercomp):
            for n in range(m, len(mclustercomp)):
                nclust = mclustercomp[n]
                nmask = np.isin(simnames, idnames[nclust])
                mmask = np.isin(simnames, idnames[mclust])
                nsimmat[m,n] = nsimmat[n,m] = np.mean(simmat[nmask][:,mmask])
        
    if '--showindividual' in sys.argv:
        profiles = sys.argv[sys.argv.index('--showindividual')+1]
        
        seqname, fsequence = masterinfo[:, masterheader.index('Motif_ID')], masterinfo[:, masterheader.index('Construct_AA_seq')]
        infname, domainseq, domainname, domainlocation = masterinfo[:, masterheader.index('Motif_ID')], masterinfo[:, masterheader.index('RBD_or_RBR_AA_Seq')], masterinfo[:, masterheader.index('Domains')], masterinfo[:, masterheader.index('Domain_Boundaries')]
        
        zscores = np.genfromtxt(profiles, dtype = str)
        znames = open(profiles, 'r').readline().split()[1:]
        kmers = zscores[:,0]
        zscores = zscores[:,1:].astype(float)
        
        keepp = []
        keepseq = []
        keepdomain = []
        for i, idname in enumerate(idnames):
            keepp.append(znames.index(idname))
            keepseq.append(list(seqname).index(idname))
            keepdomain.append(list(infname).index(idname))
        fsequence = np.array(fsequence)[keepseq]
        domainlocation = np.array(domainlocation)[keepdomain]
        znames = np.array(znames)[keepp]
        zscores = zscores[:,keepp].T
        indvnamecluster = np.arange(len(idnames),dtype = int).reshape(-1,1)
    
    #sort pwms by idnames
    spwms = []
    for i, idname in enumerate(idnames):
        spwms.append(pwms[pwmnames.index(idname)])
    # combine pwms for clusters
    if '--showindividual' in sys.argv:
        clusterpwms = spwms
    else:
        clusterpwms, clustermotifs = clustermot(spwms, idnames, mclustercomp)
    
    # assign protein names to clusters
    if '--nonames' in sys.argv:
        phylonames = []
        for i in range(len(mclustercomp)):
            phylonames.append('# '+str(i)+' ('+str(mclustlen[i])+')')
    else:
        protnames = masterinfo[:, [1,3,4]]
        if '--showindividual' in sys.argv:
            phylonames = clusterclean(idnames, protnames, indvnamecluster)
        else:
            phylonames = clusterclean(idnames, protnames, mclustercomp)
    
    # assign domaintype to protein or cluster
    dtype = []
    for i, idname in enumerate(idnames):
        di = list(masterinfo[:,1]).index(idname)
        dtype.append('-'.join(masterinfo[di,8]))
    dtype = np.array(dtype)
    
    unidtype, unidtypenum = np.unique(dtype, return_counts = True)
    domainmat = np.zeros((len(mclustercomp), len(unidtype)))
    # indicates whether different types of domains are present in motif cluster
    convergent = []
    for m, mclust in enumerate(mclustercomp):
        mdomains, mdomnum = np.unique(dtype[mclust], return_counts = True)
        domainmat[m][np.isin(unidtype, mdomains)] = mdomnum
        mdomparts = []
        for mdom in mdomains:
            mdomparts.append(np.unique(mdom.split('-')))
        dominter = reduce(np.intersect1d, mdomparts)
        if len(dominter) == 0:
            convergent.append([1])
        else:
            convergent.append([0])
    convergent = np.array(convergent)
    
    
    
    
    if '--markset' in sys.argv:
        setcutoff = float(sys.argv[sys.argv.index('--markset')+1])
        markcolor = ListedColormap(['white', 'dimgrey'])
        newset = []
        for i, idname in enumerate(idnames):
            if float(idname.strip('RNCMPT')) <= setcutoff:
                newset.append([0.])
            else:
                newset.append([1.])
        newset = np.array(newset)
        newsetcluster = []
        for m, mclust in enumerate(mclustercomp):
            newsetcluster.append([np.mean(newset[mclust])])
        newsetcluster = np.array(newsetcluster)
    
    if '--assignspecies' in sys.argv:
        Speciescolors = ListedColormap(['white', (.9,0.6,0.3,.9), (.4,0.1,0.4,.9), (.3,0.7,0.3,.9), (.6,0.2,0.2,.9),(.2,0.8,0.6,.9)])
        Acceptedspecies = ['Metazoa','Protist','Plant','Fungi','Alga']
        specieslist = masterinfo[:, [1,5]]
        specicladecolor = np.zeros((len(idnames),1))
        proteinspecies = []
        for i, idn in enumerate(idnames):
            si = list(specieslist[:,0]).index(idn)
            specicladecolor[i] = (1.+Acceptedspecies.index(specieslist[si,1]))/float(len(Acceptedspecies))
            proteinspecies.append(specieslist[si,1])
        specicladecolor = np.array(specicladecolor)
        proteinspecies = np.array(proteinspecies)
        
        specnums = np.zeros((len(mclustercomp), len(Acceptedspecies)))
        for m, mclust in enumerate(mclustercomp):
            unispec, unisnum = np.unique(proteinspecies[mclust], return_counts = True)
            for a, acs in enumerate(unispec):
                specnums[m, Acceptedspecies.index(acs)] = unisnum[a]
    else:
        Speciescolors = cm.Greys
        specnums = np.zeros((len(mclustercomp), 1))
        for m, mclust in enumerate(mclustercomp):
            specnums[m] = len(mclust)
        
    
    
    
    
    
    
    
    

    
    if '--plot' in sys.argv:
        
        if '--markset' in sys.argv:
            markgap = 0.6
        else:
            markgap = 0
        
        bottom = 0.25
        maxright = 0.7
        maxleft = 0.025
        top = 0.9
        
        pwmsize = 2
        evolsize = 0.6
        phylosize = 8
        matsize = 14
        barsize = 2
        
        mst = (maxright-maxleft)/31.
        
        
        
        fig = plt.figure(figsize = (14+0.25*len(mclustercomp[0]), len(mclustercomp)*0.3))
        ax0 = fig.add_subplot(141)
        ax0.set_position([maxleft,bottom,phylosize*mst,top-bottom])
        if '--markset' in sys.argv:
            axm = fig.add_subplot(152)
            axm.set_position([maxleft+(phylosize+pwmsize)*mst,bottom,mst*markgap,top-bottom])
        ax = fig.add_subplot(142)
        ax.set_position([maxleft+(phylosize+pwmsize+evolsize+markgap)*mst,bottom,matsize*mst,top-bottom])
        axc = fig.add_subplot(143)
        axc.set_position([maxleft+(phylosize+pwmsize+markgap)*mst,bottom,mst*evolsize,top-bottom])
        axn = fig.add_subplot(144)
        axn.set_position([maxleft+(phylosize+pwmsize+evolsize+markgap+matsize)*mst,bottom,barsize*mst,top-bottom])
        
        
        with plt.rc_context({'lines.linewidth': 3.}):
            dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = ax0)
        sorting = dn['leaves']
        ax0.tick_params(left = False, labelleft = False, right = False, labelright = False)
        ax0.spines['top'].set_visible(False)
        ax0.spines['left'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.set_xlabel(clustername)
        ax0.set_xticks(clusterticks)
        ax0.set_xticklabels(clusterticklabels)
        
        
        ax.tick_params(left = False, bottom = False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(which = 'both', left = False, bottom = False, labelleft = False)
        ax.imshow(np.sqrt(domainmat[sorting]), origin = 'lower', cmap = 'Purples', aspect = 'auto')
        ax.set_xticks(np.arange(len(unidtype)+1)-0.5, minor = True)
        ax.set_yticks( np.arange(len(phylonames)+1)-0.5, minor = True)
        ax.grid(color='silver', linewidth=0.5, which = 'minor')
        ax.set_xticks(np.arange(len(unidtype))+0.5)
        ax.set_xticklabels(unidtype, rotation = 60, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 10)
        #ax.set_yticks(np.arange(len(phylonames)))
        #ax.set_yticklabels(np.array(phylonames)[sorting]) #np.arange(len(phylonames))[sorting]) #phylonames[sorting])
        domainmat = domainmat[sorting]
        for i in range(len(domainmat)):
            for j in range(len(domainmat[0])):
                if domainmat[i, j] > 0:
                    text = ax.text(j, i, int(domainmat[i, j]),ha="center", va="center", color="white")
        
        axc.imshow(convergent[sorting],origin = 'lower', cmap = 'Reds', aspect = 'auto')
        axc.spines['top'].set_visible(False)
        axc.spines['bottom'].set_visible(False)
        axc.spines['left'].set_visible(False)
        axc.spines['right'].set_visible(False)
        axc.tick_params(which = 'both', left = False, bottom = False, labelleft = False, labelbottom = False)
        
        if '--markset' in sys.argv:
            axm.imshow(newsetcluster[sorting],origin = 'lower', cmap = 'Greys', aspect = 'auto')
            axm.spines['top'].set_visible(False)
            axm.spines['bottom'].set_visible(False)
            axm.spines['left'].set_visible(False)
            axm.spines['right'].set_visible(False)
            axm.tick_params(which = 'both', left = False, bottom = False, labelleft = False, labelbottom = False)
            for i in range(len(newsetcluster)):
                axm.text(0,i, int(np.sum(domainmat[i])*newsetcluster[sorting][i]) ,ha="center", va="center", color="white")
            axm.set_xticks([0.5])
            axm.set_xticklabels(['This study'], rotation = 60, ha = 'right', verticalalignment = 'top')
        
        for s in range(len(specnums[0])):
            Speciescolors((float(len(specnums[0]))-s)/float(len(specnums[0])))
            axn.barh(np.arange(len(specnums)), np.sum(specnums[sorting, :len(specnums[0])-s], axis = 1), color = Speciescolors((float(len(specnums[0]))-s)/float(len(specnums[0]))), edgecolor = 'k')
        axn.spines['top'].set_visible(False)
        axn.spines['left'].set_visible(False)
        axn.spines['right'].set_visible(False)
        axn.tick_params(which = 'both', left = False, labelleft = False, labelright = True)
        axn.set_ylim([-.5,len(domainmat)-.5])
        axn.set_yticks(np.arange(len(domainmat)))
        axn.set_yticklabels(phylonames[sorting], fontsize = 12)
        axn.set_xticks(np.arange(0,np.amax(np.sum(domainmat, axis = 1)),5, dtype = int))
        axn.set_xlabel('Number\nRBPs')
        # could add another bar plot on top to show how many of each type are included
        
        pheight = (0.9-bottom)/float(len(mclustercomp))
        for t, s in enumerate(sorting):
            axp = fig.add_subplot(len(mclustercomp), 1, t+1)
            axp.set_position([maxleft+phylosize*mst,bottom+pheight*(t),pwmsize*mst,pheight])
            plotpwm(clusterpwms[s], axp)
        
        
        if '--savefig' in sys.argv:
            print outname+'.jpg'
            fig.savefig(outname+'.jpg', bbox_inches = 'tight', transparent=True, dpi = dpi)
        else:
            plt.show()
        
    
    
    
    
    
    if '--showsimilarity' in sys.argv:
        strechingfactor = np.around(192./float(len(nsimmat))+ 0.5,1)
        if '--markset' in sys.argv:
            markgap = 1*strechingfactor
        else:
            markgap = 0
        
        if '--assignspecies' in sys.argv:
            specgap = 1*strechingfactor
        else:
            specgap = 0
            
        dcolors = ListedColormap(['white', 'indigo', 'firebrick', 'darkgoldenrod', 'gold', 'yellowgreen','mediumturquoise', 'seagreen', 'darkorange','sienna', 'grey'])
        recognizedtypes = ['RRM','KH', 'CCCH', 'CCHC', 'RanBP', 'CSD', 'Pumillo', 'S1', 'SAM']
        matrixaddition = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        domainmat = []
        for d, dt in enumerate(dtype):
            damt = np.zeros(6)
            for dj, dti in enumerate(dt.split('-')):
                damt[dj] = matrixaddition[recognizedtypes.index(dti)]
            domainmat.append(damt)
        domainmat = np.array(domainmat)
        if '--makelegend' in sys.argv:
            figlegend = plt.figure()
            ax = figlegend.add_subplot(111)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(which = 'both', left = False, bottom = False, labelbottom = False, labelleft = False)
            for r, rdt in enumerate(recognizedtypes):
                ax.text((r/2)*0.195+0.04, ((r+1)%2)*-0.25+0.075, rdt, fontsize=14, horizontalalignment='left', verticalalignment='center')
                ax.add_patch(Rectangle(xy=((r/2)*0.195, ((r+1)%2)*-0.25), width=0.035,height=0.15, facecolor=dcolors(matrixaddition[r]), edgecolor = 'k'))
            ax.set_xlim([-0.02,0.9])
            ax.set_ylim([0.5, -2.])
            figlegend.savefig('Domainlegend.jpg', dpi = 300, bbox_inches = 'tight')
            
            
                
        bottom = 0.1
        mst = 0.025
        pw = 1.5*strechingfactor
        phw = 3.*strechingfactor#*6
        domwidth = 3.5 *strechingfactor #*5
        matsw = 44.*mst
        matsh = 28.*mst
        maxleft = 0.025
        maxright = maxleft+(pw+phw)*mst+matsw
        top = bottom + matsh
        
        
        cmapgs = np.concatenate([[[1.,1.,1.,1.] for i in range(20)] , [cm.Oranges(0.1) for i in range(20)], [cm.Oranges(0.3) for i in range(20)], [cm.Oranges(0.6) for i in range(20)], [cm.Oranges(1.) for i in range(20)]], axis = 0)
        colmap = ListedColormap(cmapgs)
        
        
        colorbar = plt.figure(figsize = (2,0.5))
        #colmap = 'coolwarm' #'bwr'
        
        vmin = 0.
        vmax = 1.
        
        ac = colorbar.add_subplot(111)
        ac.imshow([np.linspace(vmin,vmax, 100)], origin = 'lower', cmap = colmap, aspect = 'auto', vmin = vmin, vmax = vmax)
        ac.set_xticks([0,50,99])
        
        ac.set_xticklabels([vmin, (vmin + vmax)/2., vmax])
        #ac.set_xticklabels([0,0.25,0.5,0.75,1.])
        #ac.spines['top'].set_visible(False)
        #ac.spines['left'].set_visible(False)
        #ac.spines['right'].set_visible(False)
        #ac.spines['bottom'].set_visible(False)        
        ac.tick_params(left = False, labelleft = False, right = False, labelright = False)
        #ac.set_title('Top 100 overlap')
        ac.set_title('Pearson R')
        coname = 'Pearson_colorbar.jpg' #'Top100_colorbar.jpg'
        colorbar.savefig(coname, dpi = 300,bbox_inches = 'tight')
        print 'colorbar saved'
        
        
        
        #fig = plt.figure(figsize = (len(mclustercomp)*0.07, len(mclustercomp)*0.24))
        fig = plt.figure(figsize = (len(mclustercomp)*0.07, len(mclustercomp)*0.24))
        ax0 = fig.add_subplot(121)
        ax0.set_position([maxleft,bottom,phw*mst,top-bottom])
        ax = fig.add_subplot(122)
        ax.set_position([maxleft+(pw+phw)*mst,bottom,matsw,top-bottom])
        
        with plt.rc_context({'lines.linewidth': 2.}):
            dn = dendrogram(Z, orientation = 'left', color_threshold=0, above_threshold_color='k', ax = ax0)
        sorting = dn['leaves']
        ax0.tick_params(left = False, labelleft = False, right = False, labelright = False)
        ax0.spines['top'].set_visible(False)
        ax0.spines['left'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.set_xlabel(clustername)
        ax0.set_xticks(clusterticks)
        ax0.set_xticklabels(clusterticklabels)
        
        
        ax.tick_params(left = False, bottom = False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(which = 'both', left = False, bottom = False, labelbottom = False, labelleft = False)
        im = ax.imshow(nsimmat[sorting][:, sorting][:,::-1], origin = 'lower', cmap = colmap, aspect = 'auto', vmin = vmin, vmax = vmax)
        ax.set_yticks( np.arange(len(phylonames)+1)-0.5, minor = True)
        ax.set_xticks( np.arange(len(phylonames)+1)-0.5, minor = True)
        ax.grid(color='silver', linewidth=0.5, which = 'minor')
        
        
        for s in sorting[::-1]:
            print phylonames[s], idnames[s]
        
        pheight = (top-bottom)/float(len(mclustercomp))
        for t, s in enumerate(sorting):
            axp = fig.add_subplot(len(mclustercomp), 1, t+1)
            axp.set_position([maxleft+phw*mst,bottom+pheight*(t),pw*mst,pheight])
            plotpwm(clusterpwms[s], axp)
        
        axd = fig.add_subplot(144)
        axd.set_position([maxleft+(pw+phw+0.1*domwidth)*mst+matsw,bottom,domwidth*mst,top-bottom])
        axd.spines['top'].set_visible(False)
        axd.spines['bottom'].set_visible(False)
        axd.spines['left'].set_visible(False)
        axd.spines['right'].set_visible(False)
        axd.tick_params(which = 'both', left = False, bottom = False, labelleft = False, labelbottom = False, labelright = False)
        print domainmat[sorting]
        imd = axd.imshow(domainmat[sorting], origin = 'lower', cmap = dcolors, aspect = 'auto', vmin = 0, vmax = 1)
        axd.set_yticks( np.arange(len(domainmat)+1)-0.5, minor = True)
        axd.set_xticks( np.arange(len(domainmat[0])+1)-0.5, minor = True)
        axd.grid(color='white', linewidth=1.5, which = 'minor')
        if '--markset' not in sys.argv and '--assignspecies' not in sys.argv:
            axd.tick_params(which = 'both', labelright = True)
            axd.set_yticks(np.arange(len(phylonames)))
            axd.set_yticklabels(np.array(phylonames)[sorting])
        
        if '--markset' in sys.argv:
            axm = fig.add_subplot(142)
            axm.set_position([maxleft+(pw+phw+1.1*domwidth+specgap)*mst+matsw,bottom,markgap*0.8*mst,top-bottom])
            axm.spines['top'].set_visible(False)
            axm.spines['bottom'].set_visible(False)
            axm.spines['left'].set_visible(False)
            axm.spines['right'].set_visible(False)
            axm.tick_params(which = 'both', left = False, bottom = False, labelleft = False, labelbottom = False, labelright = True)
            imm = axm.imshow(newset[sorting], origin = 'lower', cmap = markcolor, aspect = 'auto', vmin = 0, vmax = 1)
            axm.set_yticks(np.arange(len(phylonames)))
            axm.set_yticklabels(np.array(phylonames)[sorting])
            
            
        if '--assignspecies' in sys.argv:
            axs = fig.add_subplot(143)
            axs.set_position([maxleft+(pw+phw+1.1*domwidth)*mst+matsw,bottom,specgap*mst,top-bottom])
            axs.spines['top'].set_visible(False)
            axs.spines['bottom'].set_visible(False)
            axs.spines['left'].set_visible(False)
            axs.spines['right'].set_visible(False)
            axs.tick_params(which = 'both', left = False, bottom = False, labelleft = False, labelbottom = False, labelright = False)
            ims = axs.imshow(specicladecolor[sorting], origin = 'lower', cmap = Speciescolors, aspect = 'auto', vmin = 0, vmax = 1)
            if '--markset' not in sys.argv:
                axs.tick_params(which = 'both', labelright = True)
                axs.set_yticks(np.arange(len(phylonames)))
                axs.set_yticklabels(np.array(phylonames)[sorting])
        
        if '--savefig' in sys.argv:
            fig.savefig(outname+'.jpg', bbox_inches = 'tight', transparent=True, dpi = dpi)
            print 'Saved as', outname+'.jpg'
        else:
            plt.show()


 
    if '--showindividual' in sys.argv:
        
        for m, mclust in enumerate(mclustercomp):
            if len(mclust) >= 2:
                
                bottom, top, left, right = 0.1, 0.9, 0.1, 0.9
                scat = 0.7
                space = 0.1
                dom = 0.05
                ps = 0.15
                lc = len(mclust)-1
                singleh = (top-bottom)
                singlew = (right-left)
                for p in range(lc):
                    for q in range(p+1, lc+1):
                        qi = mclust[q]
                        pi = mclust[p]
                        print phylonames[pi]
                        print phylonames[qi]
                        if idmat[qi,pi] < clustercut:
                            figz = plt.figure(figsize = (6,6))
                            print (1.-idmat[qi,pi])*idmatmax, '%'
                            az = figz.add_subplot(221)
                            az.set_position([left+singlew*(ps+space),bottom+singleh*ps,singlew*scat,singleh*scat])
                            topx = np.around(topz(zscores[pi], zscores[qi]),2)
                            pearsonx = np.around(pearsonr(zscores[pi], zscores[qi])[0],2)
                            az.scatter(zscores[pi], zscores[qi], alpha = 0.5, c = 'purple', label = 'R: '+str(pearsonx)+'\n'+'Top100: '+str(topx))
                            
                            if '--equalscale' in sys.argv:
                                xlim = ylim = lim = [np.amin(zscores[[pi,qi]]), np.amax(zscores[[pi,qi]])]
                                az.set_xticks(np.arange(0,lim[-1] ,5))
                                az.set_yticks(np.arange(0,lim[-1] ,5))
                                az.set_xticklabels(np.arange(0,lim[-1] ,5, dtype = int), fontsize = 14)
                                az.set_yticklabels(np.arange(0,lim[-1] ,5, dtype = int), fontsize = 14)
                            else:
                                xlim = [np.amin(zscores[pi]), np.amax(zscores[pi])]
                                ylim = [np.amin(zscores[qi]), np.amax(zscores[qi])]
                                lim = [np.amin(zscores[[pi,qi]]), np.amin(np.amax(zscores[[pi,qi]], axis = 1))]
                                az.set_xticks(np.arange(0,xlim[-1] ,5))
                                az.set_yticks(np.arange(0,ylim[-1] ,5))
                                az.set_xticklabels(np.arange(0,xlim[-1] ,5, dtype = int), fontsize = 14)
                                az.set_yticklabels(np.arange(0,ylim[-1] ,5, dtype = int), fontsize = 14)
                            
                            
                            x100 = np.sort(zscores[pi])[-100]
                            y100 = np.sort(zscores[qi])[-100]
                            
                            az.plot(lim, lim, c = 'grey')
                            az.plot(xlim, [y100 for i in range(2)], c = 'grey', ls = '--')
                            az.plot([x100 for i in range(2)],ylim, c = 'grey', ls = '--')
                            az.spines['top'].set_visible(False)
                            az.spines['right'].set_visible(False)
                            az.set_xlabel(phylonames[pi],fontsize = 10)
                            az.set_ylabel(phylonames[qi], fontsize = 10)
                            topset = np.unique(np.append(np.argsort(zscores[pi])[-10:], np.argsort(zscores[qi])[-10:]))
                            
                            if '--annotatekmers' in sys.argv:
                                anno = '-kmers'
                                for i in topset:
                                    if zscores[pi][i] > x100 and zscores[qi][i] > y100:
                                        hl = 'left'
                                        vl = 'bottom'
                                    elif zscores[qi][i] < y100:
                                        hl = 'left'
                                        vl = 'top'
                                    else:
                                        hl = 'right'
                                        vl = 'bottom'
                                    az.annotate(kmers[i], (zscores[pi][i], zscores[qi][i]), horizontalalignment=hl, verticalalignment=vl, fontsize = 5)
                            else:
                                anno = ''
                            
                            
                            axq = figz.add_subplot(222)
                            axq.set_position([left,top-singleh*(dom+ps*0.6+space), ps*singlew, ps*singleh*0.6])
                            plotpwm(clusterpwms[qi], axq)
                            axp = figz.add_subplot(223)
                            axp.set_position([right-singlew*(dom+ps),bottom, ps*singlew, ps*singleh*0.6])
                            plotpwm(clusterpwms[pi], axp)
                            
                            
                            alignmat, alignticks, alignticklabels = align(fsequence[qi], fsequence[pi], domainlocation[qi], domainlocation[pi])
                            ay = figz.add_subplot(224)
                            ay.set_position([left+singlew*(ps+space),top-singleh*dom,singlew*scat,singleh*dom])
                            ay.imshow(alignmat, cmap = 'Greys',aspect = 'auto', vmin = 0., vmax = 1.1)
                            ay2 = ay.twiny()
                            ay2.set_position([left+singlew*(ps+space),top-singleh*dom,singlew*scat,singleh*dom])
                            ay.spines['top'].set_visible(False)
                            ay.spines['right'].set_visible(False)
                            ay.spines['bottom'].set_visible(False)
                            ay.spines['left'].set_visible(False)
                            ay2.spines['top'].set_visible(False)
                            ay2.spines['right'].set_visible(False)
                            ay2.spines['bottom'].set_visible(False)
                            ay2.spines['left'].set_visible(False)
                            ay.tick_params(left = False, labelleft = False, right = False, labelright = False, labelbottom = False, bottom = False, top = False, labeltop = True)
                            ay2.tick_params(left = False, labelleft = False, right = False, labelright = False, labelbottom = True, bottom = False, top = False, labeltop = False)
                            ay.set_xticks(alignticks[0])
                            ay.set_xticklabels(alignticklabels[0], rotation = 60, fontsize = 14)
                            ay2.set_xticks(np.arange(0, len(alignmat[0]), 100, dtype = int))
                            if not '--nolegend' in sys.argv:
                                az.legend(prop={'size': 7})
                            
                            if '--savefig' in sys.argv:
                                print 'Saved', outname+'-'+phylonames[pi].replace(' ', '').replace('(','_').replace(')','').replace('$', '').replace("\\", '')+'-'+phylonames[qi].replace(' ', '').replace('(','_').replace(')','').replace('$', '').replace("\\", '')+'-scatter_aligment'+anno+'.jpg'
                                figz.savefig(outname+'-'+phylonames[pi].replace(' ', '').replace('(','_').replace(')','').replace('$', '').replace("\\", '')+'-'+phylonames[qi].replace(' ', '').replace('(','_').replace(')','').replace('$', '').replace("\\", '')+'-scatter_aligment'+anno+'.jpg', bbox_inches = 'tight', transparent=True, dpi = dpi)
                            else:
                                plt.show()
                            plt.close()
                    



