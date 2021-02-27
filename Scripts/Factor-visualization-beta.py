# Factor-visualization.py
import numpy as np
import sys, os
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import logomaker as lm
import pandas as pd
import copy as cp
from scipy.stats import pearsonr

#### Generate pwmfigure

def ic(pw, bg):
    icout = pw*np.log2(pw/bg)
    icout[icout < 0] = 0
    icout[np.isinf(icout)] = 0
    icout[np.isnan(icout)] = 0
    return icout

def plotpfms(protpwms, nprotpwms, kpwms, nkpwms, aa, nts, combined, colorz, colorp, outname):
    ### make pandas dataframe    # 
    protframes = []
    for p, pwms in enumerate(protpwms):
        pframes = []
        for q, pwm in enumerate(pwms):
            #print pwm
            if len(pwm) != 0:
                pwmdataf = pd.DataFrame()
                for a, amino in enumerate(aa):
                    pwmdataf[amino] = pwm[:,a]
                #print pwmdataf
                pframes.append(pwmdataf)
        protframes.append(pframes)
    kmerframes = []
    for p, pwms in enumerate(kpwms):
        kframes = []
        for q, pwm in enumerate(pwms):
            if len(pwm) != 0:
                pwmdataf = pd.DataFrame()
                for a, nt in enumerate(nts):
                    pwmdataf[nt] = pwm[:,a]
                kframes.append(pwmdataf)
        kmerframes.append(kframes)
    
    lowvalpwm = False
    if len(nkpwms) > 0:
        lowvalpwm = True
        nprotframes = []
        for p, pwms in enumerate(nprotpwms):
            pframes = []
            for q, pwm in enumerate(pwms):
                #print pwm
                if len(pwm) != 0:
                    pwmdataf = pd.DataFrame()
                    for a, amino in enumerate(aa):
                        pwmdataf[amino] = pwm[:,a]
                    #print pwmdataf
                    pframes.append(pwmdataf)
            nprotframes.append(pframes)
        nkmerframes = []
        for p, pwms in enumerate(nkpwms):
            kframes = []
            for q, pwm in enumerate(pwms):
                if len(pwm) != 0:
                    pwmdataf = pd.DataFrame()
                    for a, nt in enumerate(nts):
                        pwmdataf[nt] = pwm[:,a]
                    kframes.append(pwmdataf)
            nkmerframes.append(kframes)
    
    for p, ppwms in enumerate(protframes):
        zpwms = kmerframes[p]
        nzts = len(zpwms)
        nps = len(ppwms)
        if combined:
            if nzts > 0 and nps > 0:
                fig = plt.figure('Factor'+str(p), (12, 8))
        else:
            if nps > 0:
                fig = plt.figure('Factor'+str(p)+'_p', (8, 3*nps))
            if nzts > 0:
                fig2 = plt.figure('Factor'+str(p)+'_z', (8, 3*nzts))
        if nps > 0:
            for q, pwm in enumerate(ppwms):
                #print pwm
                if combined:
                    ax = fig.add_subplot(2,nps,q+1)
                else:
                    ax = fig.add_subplot(nps,1, q+1)
                if '--infocont' in sys.argv:
                    pwm = ic(pwm, 0.05)
                lm.Logo(pwm, ax = ax, color_scheme='chemistry')
                if '--infocont' in sys.argv:
                    ax.set_yticks([0,2,4])
                    ax.set_yticklabels([0,2,4])
                    ax.set_ylim([0,np.log2(20.)])
                else:
                    ax.set_yticks([0,1])
                    ax.set_yticklabels([0,1])
                    ax.set_ylim([0,1])
                ax.tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
                ax.tick_params( axis='x',which='both', bottom=False, top=False, labelbottom=False, labeltop = False)
                for axis in ['top','bottom','left','right']:
                    if '--removespines' in sys.argv:
                        ax.spines[axis].set_visible(False)
                    else:
                        ax.spines[axis].set_linewidth(10.)
                        ax.spines[axis].set_color(colorp)
                    
                
        if nzts > 0:
            for q, pwm in enumerate(zpwms):
                if combined:
                    ax2 = fig.add_subplot(2,nzts,nzts+q+1)
                else:
                    ax2 = fig2.add_subplot(nzts,1,q+1)
                if '--infocont' in sys.argv:
                    pwm = ic(pwm, 0.25)
                lm.Logo(pwm, ax = ax2)
                if '--infocont' in sys.argv:
                    ax2.set_yticks([0,1,2])
                    ax2.set_yticklabels([0,1,2])
                    ax2.set_ylim([0,2])
                else:
                    ax2.set_yticks([0,1])
                    ax2.set_yticklabels([0,1])
                    ax2.set_ylim([0,1])
                ax2.tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
                ax2.tick_params( axis='x',which='both', bottom=False, top=False, labelbottom=False, labeltop = False)
                for axis in ['top','bottom','left','right']:
                    if '--removespines' in sys.argv:
                        ax2.spines[axis].set_visible(False)
                    else:
                        ax2.spines[axis].set_linewidth(10.)
                        ax2.spines[axis].set_color(colorz)
        if lowvalpwm:
            nppwms = nprotframes[p]
            nzpwms = nkmerframes[p]
            nnzts = len(nzpwms)
            nnps = len(nppwms)
            if combined:
                if nnzts > 0 and nnps > 0:
                    fig3 = plt.figure('Factor'+str(p)+'low', (12, 8))
            else:
                if nnps > 0:
                    fig3 = plt.figure('Factor'+str(p)+'_plow', (8, 3*nnps))
                if nnzts > 0:
                    fig4 = plt.figure('Factor'+str(p)+'_zlow', (8, 3*nnzts))
            
            if nnps > 0:
                for q, pwm in enumerate(nppwms):
                    if combined:
                        ax3 = fig3.add_subplot(2,nnps,q+1)
                    else:
                        ax3 = fig3.add_subplot(nnps,1, q+1)
                    if '--infocont' in sys.argv:
                        pwm = ic(pwm, 0.05)
                    lm.Logo(pwm, ax = ax3, color_scheme='chemistry')
                    if '--infocont' in sys.argv:
                        ax3.set_yticks([0,2,4])
                        ax3.set_yticklabels([0,2,4])
                        ax3.set_ylim([0,np.log2(20.)])
                    else:
                        ax3.set_yticks([0,1])
                        ax3.set_yticklabels([0,1])
                        ax3.set_ylim([0,1])
                    ax3.tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
                    ax3.tick_params( axis='x',which='both', bottom=False, top=False, labelbottom=False, labeltop = False)
                    for axis in ['top','bottom','left','right']:
                        if '--removespines' in sys.argv:
                            ax3.spines[axis].set_visible(False)
                        else:
                            ax3.spines[axis].set_linewidth(10.)
                            ax3.spines[axis].set_color(colorp)
            if nnzts > 0:
                for q, pwm in enumerate(nzpwms):
                    if combined:
                        ax4 = fig3.add_subplot(2,nnzts,nzts+q+1)
                    else:
                        ax4 = fig4.add_subplot(nnzts,1,q+1)
                    if '--infocont' in sys.argv:
                        pwm = ic(pwm, 0.25)
                    lm.Logo(pwm, ax = ax4)
                    if '--infocont' in sys.argv:
                        ax4.set_yticks([0,1,2])
                        ax4.set_yticklabels([0,1,2])
                        ax4.set_ylim([0,2])
                    else:
                        ax4.set_yticks([0,1])
                        ax4.set_yticklabels([0,1])
                        ax4.set_ylim([0,1])
                    ax4.tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
                    ax4.tick_params( axis='x',which='both', bottom=False, top=False, labelbottom=False, labeltop = False)
                    for axis in ['top','bottom','left','right']:
                        if '--removespines' in sys.argv:
                            ax4.spines[axis].set_visible(False)
                        else:
                            ax4.spines[axis].set_linewidth(10.)
                            ax4.spines[axis].set_color(colorz)
        #if combined:
            #fig.tight_layout()
            #fig3.tight_layout()
        if '--savefig' in sys.argv:
            print 'Saved as', outname+'_'+'Factor'+str(p)+'-motifs.jpg'
            if combined:
                if nps > 0 and nzts > 0:
                    fig.savefig(outname+'_'+'Factor'+str(p)+'-motifs.jpg', format = 'jpg', dpi = 100, bbox_inches = 'tight', transparent = True)
            else:
                if nps > 0:
                    fig.savefig(outname+'_'+'Factor'+str(p)+'-pmotifs.jpg', format = 'jpg', dpi = 100,bbox_inches = 'tight', transparent = True)
                if nzts > 0:
                    fig2.savefig(outname+'_'+'Factor'+str(p)+'-zmotifs.jpg', format = 'jpg', dpi = 100,bbox_inches = 'tight', transparent = True)
            if lowvalpwm:
                if combined:
                    if nnps > 0 and nnzts > 0:
                        fig3.savefig(outname+'_'+'Factor'+str(p)+'-lowmotifs.jpg', format = 'jpg', dpi = 100,bbox_inches = 'tight', transparent = True)
                else:
                    if nnps > 0:
                        fig3.savefig(outname+'_'+'Factor'+str(p)+'-plowmotifs.jpg', format = 'jpg', dpi = 100,bbox_inches = 'tight', transparent = True)
                    if nnzts > 0:
                        fig4.savefig(outname+'_'+'Factor'+str(p)+'-zlowmotifs.jpg', format = 'jpg', dpi = 100,bbox_inches = 'tight', transparent = True)
            #sys.exit()
        else:
            plt.show()
        
        plt.close()



#### Align sequences and generate PWMs

## Construct a PWM from pre-aligned kmers to the highest scoring k-mer
def makepwm(sed, sedscore, indexes, posins, kscore, kmotif, singles):
        fullscore = 0.
        pwwm = np.zeros((len(sed)+2*len(sed)-2,4))
        pwwmcount = np.zeros(len(sed)+2*len(sed)-2)
        for g in range(len(sed)):
            pwwm[g+len(sed)-1,singles.index(sed[g])] += sedscore
            pwwmcount[g+len(sed)-1] += 1
        # divide zscore for other positons
        pwwm[:len(sed)-1] += sedscore/4.
        pwwm[len(sed)-1+len(sed):] += sedscore/4.
        for p, i in enumerate(indexes):
            shift = posins[p]
            smot = kmotif[i]
            ssco = kscore[i]
            fullscore += ssco
            cj = 0
            for j in range(len(sed)-1+shift, len(sed)-1+shift+len(sed)):
                pwwm[j, singles.index(smot[cj])] +=ssco
                pwwmcount[j] += 1
                cj +=1
            pwwm[:len(sed)-1+shift] += ssco/4.
            pwwm[len(sed)-1+shift+len(sed):] += ssco/4.
        pwmask = pwwmcount > int(0.5*float(len(indexes)))
        pwwm = pwwm[pwmask]
        pwwm[pwwm == 0.] = 1.
        pwwm = (pwwm.T/np.sum(pwwm, axis = 1)).T
        return pwwm, fullscore

def alignshit(x, y, mcore, kmaxlength):
        pos = 0
        bestcore = 0
        bestid = 0
        for i in range(kmaxlength-mcore+1):
            spart = x[:mcore+i]
            cpart = y[kmaxlength-mcore-i:]
            ident = 0
            core = 0
            highcore = [0]
            for a in zip(spart, cpart):
                if (a[0] != '-') & (a[1] != '-'):
                    if (a[0] == a[1]) and (a[0] != 'X'):
                        ident +=1
                        core +=1
                    elif a[0] == 'X' or a[1] == 'X':
                        core +=1
                    else:
                        highcore.append(core)
                        core = 0
            if max(highcore) >= mcore or core >= mcore:
                if ident > bestid:
                    bestid = cp.copy(ident)
                    bestcore = max(max(highcore),core)
                    pos = i-kmaxlength+mcore
                elif ident == bestid and abs(i-kmaxlength+mcore) < abs(pos):
                    bestid = cp.copy(ident)
                    bestcore = max(max(highcore),core)
                    pos = i-kmaxlength+mcore
        for i in range(kmaxlength-mcore):
            spart = x[i+1:]
            cpart = y[:-1-i]
            ident = 0
            core = 0
            highcore = [0]
            for a in zip(spart, cpart):
                if (a[0] != '-') & (a[1] != '-'):
                    if (a[0] == a[1]) and (a[0] != 'X'):
                        ident +=1
                        core +=1
                    elif a[0] == 'X' or a[1] == 'X':
                        core +=1
                    else:
                        highcore.append(core)
                        core = 0
            if max(highcore) >= mcore or core >= mcore:
                if ident > bestid:
                    bestid = cp.copy(ident)
                    bestcore = max(max(highcore),core)
                    pos = i+1
                elif ident == bestid and abs(i+1) < abs(pos):
                    bestid = cp.copy(ident)
                    bestcore = max(max(highcore),core)
                    pos = i+1
        return bestid, bestcore, pos


def makeseed(nmers, nscores, mcore, kmaxlength, mident):
    seedpair = []
    seedpairscore = []
    seedpairid = []
    connectionmat = np.zeros((len(nmers), len(nmers)), dtype = np.int8)
    for mn in range(len(nmers)):
        for no in range(mn + 1, len(nmers)):
            bid, bco, bpo = alignshit(nmers[mn], nmers[no], mcore, kmaxlength)
            if bid >= mident:
                connectionmat[mn, no] = connectionmat[no, mn] = 1
            seedpairid.append(bid)
            seedpairscore.append(bid*(nscores[mn]+nscores[no]))
            seedpair.append([mn, no])
    seedpair = np.array(seedpair)
    seedpairscore = np.array(seedpairscore)
    seedpairid = np.array(seedpairid)
    if np.amax(seedpairid) >= mident:
        smask = np.where(seedpairid == np.amax(seedpairid))[0]
        spairnum = smask[np.argmax(seedpairscore[smask])]
        seedp = seedpair[spairnum]
        rest = np.copy(seedp)
        while True:
            rlen = len(rest)
            rest = np.where(np.sum(connectionmat[rest], axis = 0)>0)[0]
            if len(rest) > rlen:
                continue
            else:
                break
        rest = np.delete(rest, np.where(np.isin(rest, seedp)))
        return seedp, rest, connectionmat
    else:
        return [], [], connectionmat

def makeseedpwm(sed, sedscore, shift, kscore, kmotif, singles):
        #print sed
        #print sed, sedscore, shift, kscore, kmotif, singles
        seedscore = kscore + sedscore
        nsing = float(len(singles))
        pwwm = np.zeros((20*len(sed),len(singles)))
        pwwmcount = np.zeros(20*len(sed))
        for g in range(len(sed)):
            pwwmcount[g+10*len(sed)-len(sed)/2] += 1
            if sed[g] in singles:
                pwwm[g+10*len(sed)-len(sed)/2,singles.index(sed[g])] += sedscore
            else:
                pwwm[g+10*len(sed)-len(sed)/2,:] += sedscore/nsing
        #pwwm[:10*len(sed)-len(sed)/2] += sedscore/nsing
        #pwwm[len(sed)-10*len(sed)+len(sed)/2:] += sedscore/nsing
        cj = 0
        for j in range(10*len(sed)-len(sed)/2+shift, 10*len(sed)-len(sed)/2+shift+len(kmotif)):
            pwwmcount[j] += 1
            if kmotif[cj] in singles:
                pwwm[j,singles.index(kmotif[cj])] += kscore
            else:
                pwwm[j,:] += kscore/nsing
            cj += 1
        #pwwm[:10*len(sed)-len(sed)/2+shift] += sedscore/nsing
        #pwwm[10*len(sed)-len(sed)/2+shift+len(kmotif):] += sedscore/nsing
        return pwwm, pwwmcount, seedscore

def alignkmerpwm(apwm, acount, axscore, kmers, kmerscores, nseed, rem, conmat, singles, normalize):
    lkm = len(kmers[0])
    nsing = float(len(singles))
    #print nseed, rem
    rest = np.copy(nseed)
    #print kmers[nseed]
    while True:
        rlen = len(rest)
        orest = np.copy(rest)
        rest = np.where(np.sum(conmat[rest], axis = 0)> 0)[0]
        apwmscore = ic((apwm.T/np.sum(apwm, axis = 1)).T,0.05)
        if len(rest) > rlen:
            adkmers = rest[~np.isin(rest,orest)]
            #print kmers[adkmers]
            for ki in adkmers:
                axscore += kmerscores[ki]
                bidk = 0
                bpok = 0
                for pi in range(len(apwm)-lkm):
                    idk = 0.
                    for ri in range(lkm):
                        if kmers[ki][ri] in singles:
                            idk += apwmscore[pi+ri, singles.index(kmers[ki][ri])]
                    
                    if idk > bidk:
                        bidk = idk
                        bpok = pi
                    elif idk == bidk and abs(bpok -lkm +1) > abs(pi -lkm +1):
                        bidk = idk
                        bpok = pi
                #print bpok, bidk, kmers[ki]
                cj = 0
                for j in range(bpok, bpok + lkm ):
                    if kmers[ki][cj] in singles:
                        apwm[j,singles.index(kmers[ki][cj])] += kmerscores[ki]
                    else:
                        apwm[j,:] += kmerscores[ki]/nsing
                    acount[j] += 1
                    cj += 1
                apwmscore = ic((apwm.T/np.sum(apwm, axis = 1)).T,0.05)
                #apwm[:bpok] += kmerscores[ki]/nsing
                #apwm[bpok + lkm:] += kmerscores[ki]/nsing
        else:
            break
    #print apwm
    #print acount
    #normalize = True
    if normalize:
        # remove all positions with less than 1 count
        pwmask = acount > 0
        apwm = apwm[pwmask]
        #pseudo counts
        apwm[apwm == 0.] = .1
        #normalize by position
        apwm = (apwm.T/np.sum(apwm, axis = 1)).T
    else:
        # normalize by maximum sum, 
        apwm = apwm[acount > 0]
        #pseudo counts
        #apwm += 1.
        apwm[apwm == 0.] = .1
        apwm = apwm/np.amax(apwm) #np.sum(apwm, axis = 1))
    apwm = apwm[np.amax(apwm, axis = 1) > 0.3]
    #print apwm
    #print np.shape(apwm)
    return apwm, axscore


def cluster(kmotif, kscore, mcore, mident, singles, normalize):
    kmaxlength = len(max(kmotif, key=len))
    pwms = []
    pwmscores = []
    pwmkmersets = []
    pwmmaxscores = []
    csize = []
    bgfreq = 1./float(len(singles))
    rmsets = []
    while True:
        # determines seed pair and associated k-mers
        if len(kscore) < 2:
            break
        nseed, rem, conmat = makeseed(kmotif, kscore, mcore, kmaxlength, mident)
        # this finds best alignment of seed k-mers
        if len(nseed) < 2:
            break
        rmsets.append(kmotif[np.append(nseed,rem)])
        
        bid, bco, bpo = alignshit(kmotif[nseed[0]], kmotif[nseed[1]], mcore, kmaxlength)
        # generates pwm for seed
        print kmotif[nseed], kmotif[rem]
        seed, seedcount, seedscore = makeseedpwm(kmotif[nseed[0]], kscore[nseed[0]],bpo, kscore[nseed[1]], kmotif[nseed[1]], singles)
        if len(rem) > 0:
            # aligns rest of k-mers to the seed pwm
            pwm, nscore = alignkmerpwm(seed, seedcount, seedscore, kmotif, kscore, nseed, rem, conmat, singles, normalize)
            csize.append(len(rem)+2)
            pwmscores.append(nscore)
            pwms.append(pwm)
            pwmmaxscores.append(np.amax(kscore[np.append(nseed, rem)]))
            kmotif = np.delete(kmotif, np.append(nseed, rem), None)
            kscore = np.delete(kscore, np.append(nseed, rem), None)
            if len(kscore) < 2:
                break
        else:
            csize.append(2)
            pwmscores.append(seedscore)
            pwm = seed[seedcount > 0]
            pwm = pwm/np.sum(pwm, axis = 1)[:,None]
            pwms.append(pwm)
            pwmmaxscores.append(np.amax(kscore[nseed]))
            kmotif = np.delete(kmotif, nseed, None)
            kscore = np.delete(kscore, nseed, None)
        
    #print len(pwms), rmsets
    pwmscores = np.array(pwmscores)
    pwmmaxscores = np.array(pwmmaxscores)
    sort = np.argsort(-pwmmaxscores)
    pwmscores = pwmscores[sort]
    pwmmaxscores = pwmmaxscores[sort]
    pwms = np.array(pwms)[sort]
    csize = np.array(csize)[sort]

    return pwmscores, pwms, csize



def z_score(scores):
    divider = 1.4826*np.mean(np.absolute(scores-np.mean(scores)))
    paramet = (scores-np.mean(scores))/divider
    paramet[np.isnan(paramet)] = 0.
    paramet[np.isinf(paramet)] = 0.
    return paramet


def readin(proteinfeatures, profiles, proteinset = 'none'):
    if os.path.isfile(profiles):
        Y = np.genfromtxt(profiles)[:,1:]
        Zfeatures = np.genfromtxt(profiles, dtype = str)[:,0]
        ynames = np.array(open(profiles, 'r').readline().strip().split()[1:])
        Y = np.nan_to_num(Y).T
        pfile = True
    else:
        pfile = False
        Y = None
        ynames = None
        Zfeatures = None

    if os.path.splitext(proteinfeatures)[1] == '.npz':
        pobj = np.load(proteinfeatures)
        P = pobj['features']
        sequencefeatures = pobj['kmers']
        protnames = pobj['protnames']
        pnames = pobj['expnames']
        protnames = np.append([pnames], [protnames], axis = 0)
    else:
        pobj = open(proteinfeatures, 'r').readlines()
        protnames = []
        sequencefeatures = []
        ph = 0
        for line in pobj:
            if line[0] == '#':
                ph +=1
                protnames.append(line.strip().split()[1:])
            else:
                sequencefeatures.append(line.split()[0])
        P = np.genfromtxt(proteinfeatures, skip_header = ph)[:,1:]

    seqfeatsort = np.argsort(sequencefeatures)
    sequencefeatures = sequencefeatures[seqfeatsort]
    P = P[seqfeatsort]
    P = P.T

    if pfile:
        #sort after pnames
        if len(np.intersect1d(ynames, protnames[0])) != 0:
            protnames = np.array(protnames[0])
        else:
            protnames = np.array(protnames[1])
            
        indexes = []
        for i, ina in enumerate(protnames):
            if ina in ynames:
                indexes.append([i, list(ynames).index(ina)])
        indexes = np.array(indexes).T
        Y = Y[indexes[1]]
        P = P[indexes[0]]
        protnames = protnames[indexes[0]]
        ynames = ynames[indexes[1]]
        mprot = np.shape(P)[0]
        nprobes = np.shape(Y)[1]
    
    if os.path.isfile(proteinset):
        pset = np.genfromtxt(proteinset, dtype = str)
        pmask = np.isin(protnames, pset)
        Y = Y[pmask]
        P = P[pmask]
        protnames = protnames[pmask]
        ynames = ynames[pmask]
    
    print np.shape(Y), np.shape(P)
    return Y, ynames, Zfeatures, P, protnames, sequencefeatures






if __name__ == '__main__':
    
    profiles = sys.argv[1]
    proteinfeatures = sys.argv[2]
    proteinset =  sys.argv[3]
    maxnum = int(sys.argv[4])
    outname = os.path.splitext(profiles)[0]+'_'+os.path.splitext(os.path.split(proteinfeatures)[1])[0]

    Y, datnames, zfeatures, P, pnames, pfeatures = readin(proteinfeatures, profiles, proteinset = proteinset)
    XY = np.append(P,Y, axis = 1)
    
    if '--normP2' in sys.argv:
        P = P/np.sqrt(np.sum(P*P, axis = 1))[:, None]
    if '--normY2' in sys.argv:
        Y = Y/np.sqrt(np.sum(Y*Y, axis = 1))[:, None]
    
    # SVD for concatenated specificity and sequence k-mers
    u,s,v = np.linalg.svd(XY, full_matrices=False)
    
    # Variance explained for joint features
    s2 = s**2/np.sum(s**2)
    # eigenvector parts of joint eigenvector for z-scores
    vr = v[:, -len(Y[0]):]
    # eigenvector parts of joint eigenvector for protein features
    vp = v[:,: -len(Y[0])]

    # Rescaled explained variance captured for joint matrices
    sr2 = s**2 * np.sum(v[:, -len(Y[0]):]**2, axis = 1)/ np.sum(s**2 * np.sum(v[:, -len(Y[0]):]**2, axis = 1))
    sp2 = s**2 * np.sum(v[:,: -len(Y[0])]**2, axis = 1)/ np.sum(s**2 * np.sum(v[:, :-len(Y[0])]**2, axis = 1))
    
    # rearange vectors
    sortvec = np.argsort(-sr2)[:maxnum]
    
    u,s,v = u[:,sortvec], s[sortvec], v[sortvec]
    s2, sr2, sp2 = s2[sortvec], sr2[sortvec], sp2[sortvec]
    vr, vp = vr[sortvec], vp[sortvec]


    leny = len(zfeatures)
    lenp = len(pfeatures)
    
    colorz =[.3,1.,.7,1.]
    colorp = [.5,.8,1.,1.]
    
    
    if '--plot_scoredistribution' in sys.argv:
        for f in range(len(s)):
            fig = plt.figure(figsize = (5,2))
            ax = fig.add_subplot(111)
            #ax2 = fig.add_subplot(212)
            bins = np.linspace(np.amin(v[f]), np.amax(v[f]), 100)
            ax.hist(vr[f], bins = bins, color = colorz, alpha = 0.6, histtype = 'stepfilled', edgecolor = 'k')
            ax.hist(vp[f], bins = bins, color = colorp, alpha = 0.6, histtype = 'stepfilled', edgecolor = 'k')
            ax.set_yscale('log')
            ax.tick_params( axis='x',which='both', bottom=True, top=False, labelbottom=True, labeltop = False)
            ax.tick_params( axis='y',which='both', left=False, right=False, labelleft=False, labelright = False)
            ax.spines['left'].set_visible(False)
            ax.spines['top'].set_visible(False)
            #ax.spines['bottom'].set_visible(False)
            ax.spines['right'].set_visible(False)
     
            ax.set_xlim([np.amin(v[f]), np.amax(v[f])])

            ax.set_title('Eigenvector '+str(f+1)+' : '+str(np.around(100.*s2[f],1))+'% | R: '+str(np.around(100.*sr2[f],1))+'%| P: '+str(np.around(100.*sp2[f],1))+'%')

            plt.tight_layout()
            if '--savefig' in sys.argv:
                print 'Saved', outname+'F'+str(f)+'_score_distribution.jpg'
                fig.savefig(outname+'F'+str(f)+'_score_distribution.jpg', format = 'jpg', dpi = 200)
            else:
                plt.show()
                
            plt.close()
    
    
    
    nuctides = list(np.sort(['A', 'C', 'G', 'U']))
    aminoacid = list(np.sort(['A', 'R', 'N', 'D', 'P', 'V', 'I', 'C', 'Y', 'H', 'T', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S']))
    
    bothends = False
    if '--bothends' in sys.argv:
        bothends = True
    
    if '--makelist' in sys.argv:
        top =100
        
        names = []
        topr = []
        botr = []
        topp = []
        botp = []
        for f in range(len(s)):
            names.append(['Eigenvector'+str(f+1), 'Eigenvector'+str(f+1)+'.Z'])
            eigval = (v[f]-np.mean(v[f]))/np.std(v[f])
            eigvalr = np.around(eigval[lenp:],2)
            eigvalp = np.around(eigval[:lenp],2)
            sortr = np.argsort(eigvalr)
            sortp = np.argsort(eigvalp)
            topr.append([zfeatures[sortr[::-1][:top]],eigvalr[sortr[::-1][:top]]])
            botr.append([zfeatures[sortr[:top]],-eigvalr[sortr[:top]]])
            topp.append([pfeatures[sortp[::-1][:top]],eigvalp[sortp[::-1][:top]]])
            botp.append([pfeatures[sortp[:top]],-eigvalp[sortp[:top]]])
        np.savetxt(outname+'_eigenprofile_top100R.txt', np.append(np.concatenate(names).reshape(-1,1), np.concatenate(topr, axis = 0),axis = 1), fmt = '%s')
        np.savetxt(outname+'_eigenprofile_bot100R.txt', np.append(np.concatenate(names).reshape(-1,1),np.concatenate(botr, axis = 0),axis = 1), fmt = '%s')
        np.savetxt(outname+'_eigenprofile_top100P.txt', np.append(np.concatenate(names).reshape(-1,1),np.concatenate(topp, axis = 0),axis = 1), fmt = '%s')
        np.savetxt(outname+'_eigenprofile_bot100P.txt', np.append(np.concatenate(names).reshape(-1,1),np.concatenate(botp, axis = 0),axis = 1), fmt = '%s')
    
    # generate alignment of Z-scores:
    if '--zkmerchoice' in sys.argv:
        topn = int(sys.argv[sys.argv.index('--zkmerchoice')+1])
        minover = int(sys.argv[sys.argv.index('--zkmerchoice')+2])
        minid = int(sys.argv[sys.argv.index('--zkmerchoice')+3])
        outname += 'ztop'+str(topn)+'over'+str(minover)+'id'+str(minid)
        zkmerpwms = []
        negzkmerpwms = []
        zkmerscores = []
        zkmersize = []
        for f, factor in enumerate(v[:, lenp:]):
            print f
            #print len(factor)
            factor = factor
            factor = z_score(factor)
            sortf = np.argsort(-factor)[:topn]
            print zfeatures[sortf], factor[sortf]
            #print zfeatures[sortf]
            opwmscores, outpwms, opwmsize = cluster(zfeatures[sortf], factor[sortf], minover, minid, nuctides, True)
            #print 'Second'
            print opwmsize
            if bothends:
                sortf = np.argsort(factor)[:topn]
                print 'low'
                print zfeatures[sortf], factor[sortf]
                opwmscoresneg, outpwmsneg, opwmsizeneg = cluster(zfeatures[sortf], -factor[sortf], minover, minid, nuctides, True)
                negzkmerpwms.append(outpwmsneg)
                print opwmsizeneg
            zkmerpwms.append(outpwms)
            zkmerscores.append(opwmscores)
            zkmersize.append(opwmsize)
            
            
    if '--pkmerchoice' in sys.argv:
        topn = int(sys.argv[sys.argv.index('--pkmerchoice')+1])
        minover = int(sys.argv[sys.argv.index('--pkmerchoice')+2])
        minid = int(sys.argv[sys.argv.index('--pkmerchoice')+3])
        outname += 'ztop'+str(topn)+'over'+str(minover)+'id'+str(minid)
        pkmerpwms = []
        negpkmerpwms = []
        pkmerscores = []
        pkmersize = []
        for f, factor in enumerate(v[:, :lenp]):
            factor =factor
            factor = z_score(factor)
            sortf = np.argsort(-factor)[:topn]
            opwmscores, outpwms, opwmsize = cluster(pfeatures[sortf], factor[sortf], minover, minid, aminoacid, False)
            if bothends:
                sortf = np.argsort(factor)[:topn]
                opwmscoresneg, outpwmsneg, opwmsizeneg = cluster(pfeatures[sortf], -factor[sortf], minover, minid, aminoacid, False)
                negpkmerpwms.append(outpwmsneg)
            pkmerpwms.append(outpwms)
            pkmerscores.append(opwmscores)
            pkmersize.append(opwmsize)
    
    
    if '--infocont' in sys.argv:
        outname += '-infocont'
    
    combined = False
    if '--combined' in sys.argv:
        combined = True
    

    plotpfms(pkmerpwms, negpkmerpwms, zkmerpwms, negzkmerpwms, aminoacid, nuctides, combined, colorz, colorp, outname)







