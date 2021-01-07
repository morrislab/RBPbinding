# program to generate PFMs from RNAcompete zscore profiles
# run with
# python ~/Work/RBPs/Scripts/Motif_finder/extract_motifs_zscores.py --zscores ../Zscores_426_origrncmpt.txt --outname Zscores_426_origrncmpt_pwm_context_top100 allpwm --kmerchoice top 100 --clusterparams 5 4 4 2

import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from operator import itemgetter
import operator
from Bio import pairwise2
import copy as cp 
from sklearn.externals import joblib
from scipy import stats




class readinparms:        
    def __init__(self):
        self.inputparams = {}
        self.inputvecs = []
        
        ## provide 7-mer RNA sequence specificity
        if '--zscores' in sys.argv: # measured Zscores
            multiple = True
            facfile = sys.argv[sys.argv.index('--zscores')+1]
            facend = os.path.splitext(facfile)[-1]
            if facend == '.txt':
                paramet = np.genfromtxt(facfile, dtype = str)
                parnames = paramet[:,0]
                paramet = paramet[:,1:].astype(float)
                protnames = np.array(open(facfile, 'r').readline().strip().split()[1:])
                #self.paramet = paramet
            elif facend == '.npy' or facend == '.npz':
                npload = np.load(facfile)
                paramet = npload['motifs']
                parshape = np.shape(paramet)
                if parshape[0] < parshape[1]:
                    paramet = paramet.T
                    print np.shape(paramet)
                parnames = npload['kmers']
                protnames = npload['protnames']
            self.inputparams['protnames'] = protnames
        elif '--predictedfile' in sys.argv: # predicted zscores
            multiple = True
            facfile = sys.argv[sys.argv.index('--predictedfile')+1]
            facend = os.path.splitext(facfile)[-1]
            if facend == '.txt':
                paramet = np.genfromtxt(facfile, dtype = str)
                parnames = paramet[:,0]
                paramet = paramet[:,1:].astype(float)
                protnames = np.array(open(facfile, 'r').readline().strip().split()[1:])
                #self.paramet = paramet
            elif facend == '.npy' or facend == '.npz':
                npload = np.load(facfile)
                #print npload.files
                paramet = npload['profiles']
                parshape = np.shape(paramet)
                if parshape[0] < parshape[1]:
                    paramet = paramet.T
                    print np.shape(paramet)
                parnames = npload['kmers']
                protnames = npload['names']
            if '--cutpname' in sys.argv:
                nprname = []
                for prot in protnames:
                    nprname.append(prot.split('|')[0])
                protnames = np.array(nprname)
            self.inputparams['protnames'] = protnames.astype(str)
        elif '--predictlight' in sys.argv: # top100 file for predicted zscores
            multiple = True
            facfile = sys.argv[sys.argv.index('--predictlight')+1]
            paramet = np.genfromtxt(facfile, dtype = str)
            parnames = paramet[:,np.arange(0, len(paramet[0]), 2 ,dtype = int)]
            parnames = np.unique(np.concatenate(parnames))
            protnames = np.array(open(facfile, 'r').readline().strip().split()[1:])[np.arange(0, len(paramet[0]), 2 ,dtype = int)]
            #print protnames
            paramat = np.zeros((len(parnames), len(protnames)))
            for p, prot in enumerate(protnames):
                sortk = np.argsort(paramet[:,p*2])
                #print paramet[:,p*2][sortk], paramet[:,p*2+1][sortk]
                paramat[np.isin(parnames, paramet[:,p*2][sortk]), p] = paramet[sortk,p*2+1].astype(float)
            paramet = paramat
            self.inputparams['protnames'] = protnames.astype(str)    
        else:
            print "Please provide a file"
            sys.exit()
       
        # z-score profiles if not automatically done
        if '--normalize_zscores' in sys.argv:
            # robust variance for normalization
            divider = 1.4826*np.mean(np.absolute(paramet-np.mean(paramet, axis =0)),axis = 0)
            divider[divider == 0] = 1
            paramet = (paramet-np.mean(paramet, axis =0))/divider
        
        self.inputvecs.append(paramet)
        self.inputvecs.append(np.array(parnames, dtype = str))
        self.inputparams['add'] = multiple
        
        
        if '--outname' in sys.argv:
            outname = sys.argv[sys.argv.index('--outname')+1] # name of file
            outype = sys.argv[sys.argv.index('--outname')+2] # highpwm or allpwm defines if secondary motifs are included
            self.inputparams['outname'] = outname
            self.inputparams['outtype'] = outype
        else:
            print "Please provide outname"
            sys.exit()
            
        if '--kmerchoice' in sys.argv:
            kcuttype = sys.argv[sys.argv.index('--kmerchoice')+1] # top
            kmercut = float(sys.argv[sys.argv.index('--kmerchoice')+2]) # 10 for example
            self.inputparams['kmercut'] = kmercut
            self.inputparams['kcuttype'] = kcuttype
        else:
            print 'Please define how to choose kmer set'
            sys.exit()
        
        if '--clusterparams' in sys.argv:
            cluster = True
            minident = int(sys.argv[sys.argv.index('--clusterparams')+1]) # minimum identity that 7-mers nedd to be in one pwm
            mincore = int(sys.argv[sys.argv.index('--clusterparams')+2]) # minimum overlap of 7-mers in alignment
            minstretch = int(sys.argv[sys.argv.index('--clusterparams')+3]) # minimum continous identical stretch of two k-mers to be aligned
            minkmers = int(sys.argv[sys.argv.index('--clusterparams')+4]) # minumum k-mers for position in PFM to be included
            self.inputparams['minident'] = minident
            self.inputparams['mincore'] = mincore
            self.inputparams['minstretch'] = minstretch
            self.inputparams['minkmer'] = minkmers
        else:
            cluster = False
        self.inputparams['cluster'] = cluster    
        
        

def kmerset(cuttype, cutoff, kvals, ktxt):
#### in case kmers are all same length
    if cuttype == "top%":
        percint = int(float(len(kvals))*float(cutoff))
        mask = np.argsort(-kvals)
        kmerscores = kvals[mask][:percint]
        kmermotifs = ktxt[mask][:percint]
        kmermotifs = np.array(kmermotifs, dtype = 'object' )
        kmerscores = np.array(kmerscores)
        
    elif cuttype == "top":
        topint = int(cutoff)
        mask = np.argsort(-kvals)
        kmerscores = kvals[mask][:topint]
        kmermotifs = ktxt[mask][:topint]
        kmermotifs = np.array(kmermotifs, dtype = 'object' )
        kmerscores = np.array(kmerscores)

    elif cuttype == "threshold":
        thresh = float(cutoff)
        mask = kvals>thresh
        kmerscores = kvals[mask]
        kmermotifs = ktxt[mask]
        kmermotifs = np.array(kmermotifs, dtype = 'object' )
        kmerscores = np.array(kmerscores)
        
    elif cuttype == "thresholdmax":
        thresh = max(kvals)*float(cutoff)
        mask = kvals>thresh
        kmerscores = kvals[mask]
        kmermotifs = ktxt[mask]
        kmermotifs = np.array(kmermotifs, dtype = 'object' )
        kmerscores = np.array(kmerscores)
    else:
        print 'kmerchoice not understood!'
        sys.exit()
    #print kmermotifs
    return kmerscores, kmermotifs


def alignkmer(x, y, mcore, minstretch, kmaxlength):

    pos = 0
    bestcore = 0
    bestid = 0
    for i in range(kmaxlength-mcore+1):
        spart = x[:mcore+i]
        cpart = y[kmaxlength-mcore-i:]
        cmat = [a[0] == a[1] for a in zip(spart, cpart)]
        ident = 0
        core = 0
        highcore = [0]
        for cma in cmat:
            if cma == True:
                ident +=1
                core +=1
            elif cma == False:
                highcore.append(core)
                core = 0
        if max(highcore) >= minstretch or core >= minstretch:
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
        cmat = [a[0] == a[1] for a in zip(spart, cpart)]
        ident = 0
        core = 0
        highcore = [0]
        for cma in cmat:
            if cma == True:
                ident +=1
                core +=1
            elif cma == False:
                highcore.append(core)
                core = 0
        if max(highcore) >= minstretch or core >= minstretch:
            if ident > bestid:
                bestid = cp.copy(ident)
                bestcore = max(max(highcore),core)
                pos = i+1
            elif ident == bestid and abs(i+1) < abs(pos):
                bestid = cp.copy(ident)
                bestcore = max(max(highcore),core)
                pos = i+1
    
    return bestid, bestcore, pos


def makeseed(nmers, nscores, mcore, minstretch, kmaxlength, mident):
    seedpair = []
    seedpairscore = []
    seedpairid = []
    pairidmat = np.ones((len(nmers),len(nmers)))*kmaxlength
    for mn in range(len(nmers)):
        for no in range(mn + 1, len(nmers)):
            bid, bco, bpo = alignkmer(nmers[mn], nmers[no], mcore, minstretch, kmaxlength)
            pairidmat[mn,no] = pairidmat[no,mn] = bid
            seedpairid.append(bid)
            seedpairscore.append(bid*(nscores[mn]+nscores[no]))
            seedpair.append([mn, no])
    seedpair = np.array(seedpair)
    seedpairscore = np.array(seedpairscore)
    seedpairid = np.array(seedpairid)
    smask = np.where(seedpairid == np.amax(seedpairid))[0]
    spairnum = smask[np.argmax(seedpairscore[smask])]
    seedp = seedpair[spairnum]
    rest = np.copy(seedp)
    lnow = 0
    while True:
        rest = np.append(rest, np.where(np.amax(pairidmat[rest],axis = 0) >= mident)[0])
        rest = np.unique(rest)
        if len(rest) == len(nmers) or len(rest) == lnow:
            break
        lnow = len(rest)
    leftover = np.delete(np.arange(len(nmers), dtype = int), rest)
    rest = np.delete(rest, seedp).astype(int)
    return seedp, rest, leftover

def makeseedpwm(sed, sedscore, shift, kscore, kmotif):
        print print_alignment(sed, 6), sedscore
        print print_alignment(kmotif, 6+shift), kscore
        seedscore = kscore + sedscore
        pwwm = np.zeros((len(sed)+2*len(sed)-2,4))
        pwwmcount = np.zeros(len(sed)+2*len(sed)-2)
        for g in range(len(sed)):
            pwwm[g+len(sed)-1,singles.index(sed[g])] += sedscore
            pwwmcount[g+len(sed)-1] += 1
        pwwm[:len(sed)-1] += sedscore/4.
        pwwm[len(sed)-1+len(sed):] += sedscore/4.
        cj = 0
        for j in range(len(kmotif)-1+shift, len(kmotif)-1+shift+len(kmotif)):
            pwwm[j,singles.index(kmotif[cj])] += kscore
            pwwmcount[j] += 1
            cj += 1
        pwwm[:len(kmotif)-1+shift] += sedscore/4.
        pwwm[len(kmotif)-1+shift+len(kmotif):] += sedscore/4.
        return pwwm, pwwmcount, seedscore

def print_alignment(ksds, offs):
    offleft = ''
    offright = ''
    for o in range(offs):
        offleft+= '-'
    for o in range(len(ksds)+offs,len(ksds)*3-2):
        offright +='-'
    return offleft+ksds+offright
    


def alignkmerpwm(apwm, acount, axscore, kmers, kmerscores,minkmer):
    lkm = len(kmers[0])
    for ki in range(len(kmers)):
        axscore += kmerscores[ki]
        bidk = 0
        bpok = 0
        for pi in range(len(apwm)-lkm):
            idk = 0.
            for ri in range(lkm):
                idk += apwm[pi+ri, singles.index(kmers[ki][ri])]
            if idk > bidk:
                bidk = idk
                bpok = pi
            elif idk == bidk and abs(bpok -lkm +1) > abs(pi -lkm +1):
                bidk = idk
                bpok = pi
        print print_alignment(kmers[ki], bpok), kmerscores[ki]
        cj = 0
        for j in range(bpok, bpok + lkm ):
            apwm[j,singles.index(kmers[ki][cj])] += kmerscores[ki]
            acount[j] += 1
            cj += 1
        apwm[:bpok] += kmerscores[ki]/4.
        apwm[bpok + lkm:] += kmerscores[ki]/4.
    pwmask = acount >= minkmer
    apwm = apwm[pwmask]
    apwm[apwm == 0.] = 1.
    #print apwm
    apwm = (apwm.T/np.sum(apwm, axis = 1)).T
    return apwm, axscore


def cluster(kmotif, kscore, mcore, mident, minstretch, minkmer):
    kmaxlength = len(max(kmotif, key=len))
    pwms = []
    pwmscores = []
    csize = []
    while True:
        # find seed as best matching 7-mers with highest z-score, determine which other kmers need to be included into current pwm and which wont
        nseed, rem, lefto = makeseed(kmotif, kscore, mcore, minstretch, kmaxlength, mident)
        #print nseed, rem, lefto, mident
        if len(rem) >= minkmer -2 and len(rem) > 0:
            print 'Align', len(rem)+2
            csize.append(len(rem)+2)
            # determine best alignment of seed 7-mers
            bid, bco, bpo = alignkmer(kmotif[nseed[0]], kmotif[nseed[1]], mcore, minstretch, kmaxlength)
            # generate seed pwm from seed 7-mers
            seed, seedcount, seedscore = makeseedpwm(kmotif[nseed[0]], kscore[nseed[0]],bpo, kscore[nseed[1]], kmotif[nseed[1]])
            # add other 7-mers to seed pwm
            pwm, nscore = alignkmerpwm(seed, seedcount, seedscore, kmotif[rem], kscore[rem], minkmer)

            pwmscores.append(nscore)
            pwms.append(pwm)
            if len(lefto) < minkmer:
                break
            else:
                kmotif = kmotif[lefto]
                kscore = kscore[lefto]
        else:
            break
    
    pwmscores = np.array(pwmscores)
    sort = np.argsort(-pwmscores)
    pwmscores = pwmscores[sort]
    pwms = np.array(pwms)[sort]
    csize = np.array(csize)[sort]

    return pwmscores, pwms, csize

class writeout:
    def __init__(self,writeouttype, outname):
        print writeouttype
        if writeouttype == 'allPWM' or writeouttype == 'allpwm':
            self.objectname = os.path.splitext(outname)[0]+'.almot' 
            self.wobj = open(self.objectname, 'w')
            
        elif writeouttype == 'highPWM' or writeouttype == 'highpwm':
            self.objectname = os.path.splitext(outname)[0]+'.hmot'
            self.wobj = open(self.objectname, 'w')
        self.wobj.close()
        print self.objectname, 'cleared'
    
    def addmotif(self,writeouttype, pwms, pwmscores, csize, hname):
        wobj = open(self.objectname, 'a+')
        wobj.write('Motif\t'+ hname+addw+'\n')
        if '--noprint' not in sys.argv:
            print hname
        if writeouttype == 'allPWM' or writeouttype == 'allpwm':
            for p, pwm in enumerate(pwms):
                main_logo = ''
                for line in pwm:
                    if np.any(line > 0.5):
                        main_logo += singles[np.argmax(line)]
                    else:
                        main_logo += 'N'
                if '--noprint' not in sys.argv:
                    print '#\t', main_logo, csize[p], pwmscores[p]
                    print 'Pos\tA\tC\tG\tU'
                wobj.write('#'+main_logo+' '+str(csize[p])+' '+str(pwmscores[p])+'\n'+'Pos A C G U\n')
                for i in range(len(pwm)):
                    if '--noprint' not in sys.argv:
                        print i+1, ' '.join(np.around(pwm[i],6).astype(str))
                    wobj.write(str(i+1)+' '+' '.join(np.around(pwm[i],6).astype(str))+'\n')
            
        if writeouttype == 'highPWM' or writeouttype == 'highpwm':
            #print pwms
            hpwm = np.around(pwms[0], 6)
            #print hpwm
            main_logo = ''
            for line in hpwm:
                if np.any(line > 0.5):
                    main_logo += singles[np.argmax(line)]
                else:
                    main_logo += 'N'
            if '--noprint' not in sys.argv:
                print '#\t', main_logo, csize[0], pwmscores[0]
                print 'Pos\tA\tC\tG\tU'
            wobj.write('#\t'+main_logo+' '+str(csize[0])+' '+str(pwmscores[0])+'\n'+'Pos\tA\tC\tG\tU\n')
            for i in range(len(hpwm)):
                if '--noprint' not in sys.argv:
                    print i+1, '\t'.join(hpwm[i].astype(str))
                wobj.write(str(i+1)+'\t'+'\t'.join(np.around(hpwm[i],6).astype(str))+'\t\n')
        wobj.write('\n\n')
        wobj.close()
        
            #### Nucleotide abbreviations:
            # R: A, G
            # Y: C, U
            #W 	Weak (A or T)
            # S 	Strong (G or C)
            #M 	Amino (A or C)
            #K 	Keto (G or T)
            #B 	Not A (G or C or T)
            #H 	Not G (A or C or T)
            #D 	Not C (A or G or T)
            #V 	Not T (A or G or C)




if __name__ == '__main__':
    singles = ['A', 'C', 'G', 'U'] 
    if '--dna' in sys.argv:
        singles = ['A', 'C', 'G', 'T']
    readin = readinparms()
    if '--addword' in sys.argv:
        addw = sys.argv[sys.argv.index('--addword')+1]
    else:
        addw = ''
    
    parameter = readin.inputparams
    writeinto = writeout(parameter['outtype'], parameter['outname'])
    #ppar = readin.inputvecs[0]
    
    if parameter['add']:
        ffnum = len(readin.inputvecs[0][0])
        
        for f in range(ffnum):
            #if parameter['protnames'][f] == 'RNCMPT00272':
                kscore, kmot = kmerset(parameter['kcuttype'], parameter['kmercut'], readin.inputvecs[0][:,f], readin.inputvecs[1])
                #print kscore, kmot
                pscores, ps, cs = cluster(kmot, kscore, parameter['mincore'], parameter['minident'], parameter['minstretch'], parameter['minkmer'])
                #print pscores, ps, cs
                writeinto.addmotif(parameter['outtype'],ps, pscores, cs, parameter['protnames'][f])
    else:
            #print 'ppar is', np.argsort(readin.inputvecs[0][:,f])
            kscore, kmot = kmerset(parameter['kcuttype'], parameter['kmercut'], readin.inputvecs[0], readin.inputvecs[1])
            pscores, ps, cs = cluster(kmot, kscore, parameter['mincore'], parameter['minident'])
            #print pscores, ps, cs
            writeinto.addmotif(parameter['outtype'],ps, pscores, cs, parameter['protnames'])
            


    
