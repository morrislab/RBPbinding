import numpy as np
import sys, os 
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

nucleotideabbreviations = np.array(['A','C','G','U','M','R','W','S','Y','K','V','H','D','B','N'])
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
scormotmat = np.copy(nucabbmat)
scormotmat[-1,:] = 1.
scormotmat = (scormotmat.T/np.sum(scormotmat, axis = 1)).T

def msim(word1, word2):
    prodsim = 0.
    for le in range(len(word1)):
        prodsim += np.sum(scormotmat[nucleotideabbreviations == word1[le]]* scormotmat[nucleotideabbreviations == word2[le]])
    return prodsim

# measured similarity of IUPAC motifs
def motifsimilarity(mot1, motset2):
    sim = 0.
    alseq = []
    mout = -1
    for m, mot2 in enumerate(motset2):
        if len(mot1) < len(mot2):
            motA = mot1
            motB = mot2
        else:
            motA = mot2
            motB = mot1
        for g in range(len(motB)-len(motA)+1):
            tsim = msim(motB[g:g+len(motA)], motA)
            #print motB[g:g+len(motA)], motA, tsim
            if tsim > sim:
                sim = np.copy(tsim)
                alseq = [motB[g:g+len(motA)], motA]
                mout = m
        for g in range(1, len(motA)-2):
            ma = motA[g:]
            mb = motB[:len(ma)]
            tsim = msim(mb, ma)
            #print ma, mb, tsim
            if tsim > sim :
                sim = np.copy(tsim)
                alseq = [mb, ma]
                mout = m
            ma = motA[:-g]
            mb = motB[-len(ma):]
            tsim = msim(mb, ma)
            #print ma, mb, tsim
            if tsim > sim:
                sim = np.copy(tsim)
                alseq = [ma, mb]
                mout = m
    return sim, alseq, mout


def readmat(matfile):
    names = []
    sequences = []
    interface = []
    distance = []
    obj = open(matfile, 'r').readlines()
    for l, line in enumerate(obj):
        if '>' == line[0]:
            names.append(line[1:].strip())
            sequences.append(obj[l+1].strip().replace(',',''))
            interface.append(np.array(obj[l+2].strip().split(','), dtype = int))
            distance.append(np.array(obj[l+6].strip().split(','), dtype = str))
    return np.array(names), np.array(sequences), np.array(interface), np.array(distance)

def readstructure(fasta):
    fasta = open(fasta,'r').readlines()
    names = []
    sequences = []
    seconstr = []
    collectseq = False
    collectstr = False
    for l, line in enumerate(fasta):
        if line[0] == '>':
            line = line[1:].strip().split(':')
            if collectseq:
                sequences.append(sequence)
            if collectstr:
                seconstr.append(struc.replace(' ', 'C'))
            if line[-1] == 'sequence':
                names.append(line[0].lower()+'_'+line[1])
                collectseq = True
                collectstr = False
                sequence = ''
            elif line[-1] == 'secstr':
                collectseq = False
                collectstr = True
                struc = ''
                
        elif collectseq:
            sequence += line.strip('\n')
        elif collectstr:
            struc += line.strip('\n')
    if collectseq:
            sequences.append(sequence)
    if collectstr:
            seconstr.append(struc.replace(' ', 'C'))
    return np.array(names), np.array(sequences), np.array(seconstr)
            
def readfasta(fasta):
    fasta = open(fasta,'r').readlines()
    names = []
    sequences = []
    for l, line in enumerate(fasta):
        if line[0] == '>':
            names.append(line[1:].strip())
            sequences.append(fasta[l+1].strip())
    return np.array(names), np.array(sequences)
            
def readinmotif(tmot):
    obj = open(tmot, 'r').readlines()
    names = []
    iupac = []
    chains = []
    for l, line in enumerate(obj):
        names.append(line.strip().split()[0][1:])
        cs = 0
        iup = []
        chain = []
        for sl, s in enumerate(line):
            if s == '#':
                cs +=1
            if s != '#' and cs > 0:
                if cs == 1:
                    piup = line[sl:].split(',##')[0].split(',')
                    iup.append(piup[:len(piup)/2])
                if cs ==2:
                    chain.append(line[sl+1])
                cs = 0
        ochain = []
        for p, piup in enumerate(iup):
            ochain.append([chain[p] for s in range(len(piup))])
        iupac.append(np.concatenate(iup))
        chains.append(np.concatenate(ochain))

    return np.array(names), np.array(iupac), np.array(chains)

def removeends(al1, al2):
    start = 0
    end = len(al1)
    mask = (np.array(list(al1))!='-') * (np.array(list(al2))!='-')
    check = np.cumsum(mask)-np.cumsum(~mask)
    ed = False
    for i in range(len(al1)):
        if check[i] < 0:
            ed = True
        if ed and mask[i]:
            start = i
            break
    check = np.cumsum(mask[::-1])-np.cumsum(~mask[::-1])
    ed = False
    for i in range(len(al1)):
        if check[i] <0:
            ed = True
        if ed and mask[-i-1]:
            end = len(al1) - i
            break
    return start, end
    
    
def align(al1, al2, vec2, delim, cutends = None):
    # remove ends that have only a few residues followed by long stretches of ----
    #print al1 
    #print al2
    #print vec2
    if cutends is not None:
        st, en = removeends(al1, al2)
    else:
        st, en = 0, len(al1)
    start = np.sum(np.array(list(al1))[:st]!='-')
    j = np.sum(np.array(list(al2))[:st]!='-')
    al1 = al1[st:en]
    al2 = al2[st:en]
    
    # compute local identity dessity in window of 12
    if cutends is not None:
        localid = np.zeros(len(al1))
        for i, al in enumerate(al2):
            if al == al1[i]:
                localid[i] = 1.
        locid = np.copy(localid)
        for i in range(len(localid)):
            locid[i]=max(np.sum(localid[i:min(len(locid),i+10)])/10.,np.sum(localid[max(0,i-10):i])/10.)
        
        # get new start and end points for identical regions, remove flanking regions with low identity
        locid =np.where(locid >= cutends)[0]
        if len(locid) > 0:
            nst = locid[0]
            nen = locid[-1]
        else:
            return 'X', 'X', 0., [0,0]
    else:
        nst = 0
        nen = len(al1)
    # determine the start and end location in the sequence of al1 as reference
    start += np.sum(np.array(list(al1))[:nst]!='-')
    end = start + np.sum(np.array(list(al1))[:nen]!='-')
    j += np.sum(np.array(list(al2))[:nst]!='-')
    
    al1 = al1[nst:nen]
    al2 = al2[nst:nen]
    
    #print al1
    #print al2
    #print j, vec2[j]
    ## align the given score of the seconda alignment to the sequence of the first
    cover = 0
    out = ''
    outseq = ''
    for i, al in enumerate(al2):
        if al != '-':
            if al1[i] != '-':
                out += str(vec2[j])+delim
                outseq += al
                if al == al1[i]:
                    cover += 1.
            j += 1
    return out.strip(delim), outseq, cover, [int(start), int(end)]





  
def combine_pdb(motifcomb, chaincomb, seqcomb, facecomb, distancecomb, namescomb):
    lencomb = []
    for n, nc in enumerate(namescomb):
        lencomb.append(len(seqcomb[n]))
        chaincomb[n] = chaincomb[n].astype('|S15')
        for m, mot in enumerate(chaincomb[n]):
            chaincomb[n][m] = nc+"."+mot
    motifcomb = np.concatenate(motifcomb)
    chaincomb = np.concatenate(chaincomb)
    seed = np.argmax(lencomb)
    namescomb = namescomb[seed]
    seedseq = seqcomb[seed]
    outfacecomb = np.zeros((len(seqcomb), len(seedseq)))
    outdistancecomb = np.zeros((len(seqcomb), len(seedseq)))
    for s, seq in enumerate(seqcomb):
        alignment = pairwise2.align.globalds(seedseq, seq, matrix, -11,-1)
        mask = np.array(list(alignment[0][0])) != '-'
        mask2 = np.array(list(alignment[0][1])) != '-'
        cface = np.zeros(len(alignment[0][0]), dtype = int)
        cface[mask2] = facecomb[s]
        outfacecomb[s] = cface[mask]
        cface = np.ones(len(alignment[0][0]))*100.
        cface[mask2] = distancecomb[s]
        outdistancecomb[s] = cface[mask]
    outdistancecomb = np.amin(outdistancecomb, axis = 0)
    outfacecomb = np.amax(outfacecomb, axis = 0)
    
    return motifcomb, chaincomb, seedseq, outfacecomb, outdistancecomb, namescomb
    
def coverage(seq1, seq2):
    c = 0
    for s in range(len(seq1)):
        if seq1[s] == seq2[s]:
            c += 1
    return c
    
            

if __name__ == '__main__':

    idfile = np.load(sys.argv[1])
    idcut = float(sys.argv[2])

    tempname, tempmotif, tempchain = readinmotif(sys.argv[3])
    motifile = np.genfromtxt(sys.argv[4], dtype = str)
    motsimcut = float(sys.argv[5])
    
    rncmptlist = np.genfromtxt(sys.argv[6], dtype = str)
    
    intnames, intseq, interface, intdistance = readmat(sys.argv[7])
    fastaname, fastaseq = readfasta(sys.argv[8])
    secname, secseq, seqstruc = readstructure(sys.argv[9])
    
    pdbidfile = np.load(sys.argv[10])
    proteinnames = np.genfromtxt(sys.argv[11], dtype = str)
    
    outname = sys.argv[12]
    
    
    outname = outname+'-id'+str(idcut)+'-msim'+str(motsimcut)+'.txt'

    
    pidents = pdbidfile['identmat']
    pdbnames = pdbidfile['names']
    
    
    # align identity matrix to interfaces and 
    keep = []
    sort = []
    sortt = []
    for p, pdb in enumerate(pdbnames):
        if '__' in pdb:
            pdbnames[p] = pdb.split('__')[0]
        if pdbnames[p] in intnames and pdbnames[p] in tempname:
            sort.append(list(intnames).index(pdbnames[p]))
            sortt.append(list(tempname).index(pdbnames[p]))
            keep.append(p)
    intnames = intnames[sort]
    intseq = intseq[sort]
    interface = interface[sort]
    intdistance = intdistance[sort]
    tempmotif = tempmotif[sortt]
    tempname = tempname[sortt]
    tempchain = tempchain[sortt]
    pidents = pidents[keep][:,keep]
    pdbnames = pdbnames[keep]
    
    # generate list of unique pdbs and determine common interface mask by averaging over all masks (take min dist) and use longest sequence on both ends
    used = []
    # make unique pdb masks
    # if two pds have more than 70% identity, align the sequences and take the closest distance of the aligned residues as interface information
    intnames1 = []
    intseq1 = []
    interface1 = []
    intdistance1 = []
    tempmotif1 = []
    tempname1 = []
    tempchain1 = []
    print len(pdbnames)
    for p, pdb in enumerate(pdbnames):
        if p not in used:
            combine = np.where(pidents[p] >= 70.)[0]
            print "Combine", pdb, pdbnames[combine]
            #print -np.sort(-pidents[p])[:10], pdbnames[np.argsort(-pidents[p])[:10]]
            if len(combine) > 1:
                nmotif, nchain, nseq, nint, ndist, nname = combine_pdb(tempmotif[combine], tempchain[combine], intseq[combine], interface[combine], intdistance[combine], intnames[combine])
                
            else:
                nmotif, nchain, nseq, nint, ndist, nname = tempmotif[p], tempchain[p], intseq[p], interface[p], intdistance[p], intnames[p]
            intnames1.append(nname)
            intseq1.append(nseq)
            interface1.append(nint)
            intdistance1.append(ndist)
            tempmotif1.append(nmotif)
            tempname1.append(nname)
            tempchain1.append(nchain)
            used = np.append(used, combine)
    print 'Left', len(intnames1)
    print '\n'.join(np.sort(intnames1))
    
    
    intnames,intseq,interface,intdistance,tempmotif,tempname,tempchain = np.array(intnames1), np.array(intseq1),np.array(interface1),np.array(intdistance1),np.array(tempmotif1),np.array(tempname1),np.array(tempchain1)
    tempname = np.copy(intnames)
    
    # load identities to RNCMPT sequences
    idents = idfile['identmat']
    pdbnames = idfile['names2']
    rrmnames = idfile['names']
    
    # sort all RNCMPT features to identity matrix
    proteinidname = []
    sort_seq = []
    for r, rrmname in enumerate(rrmnames):
        proteinidname.append(rrmname.split('__')[0].split('||')[1])
        if '>' in rrmname:
            rrmnames[r] = rrmname[1:]
        sort_seq.append(list(fastaname).index(rrmnames[r]))
    proteinidname = np.array(proteinidname)
    mask = np.isin(proteinidname,rncmptlist)
    rrmnames = rrmnames[mask]
    proteinidname = proteinidname[mask]
    idents = idents[:,mask]
    fastaname = fastaname[sort_seq][mask]
    fastaseq = fastaseq[sort_seq][mask]
    
    
    
    for p, pdbname in enumerate(pdbnames):
        if '>' in pdbname:
            pdbnames[p] = pdbname[1:]
        if '__' in pdbname:
            pdbnames[p] = pdbname.split('__')[0]
    
    sort_pdb = []
    for p, pdbname in enumerate(intnames):
        sort_pdb.append(list(pdbnames).index(pdbname))
    pdbnames = pdbnames[sort_pdb]
    idents = idents[sort_pdb]
     
    outobj = open(outname,'w')
    outobj2 = open(os.path.splitext(outname)[0]+'_all.txt', 'w')
    
    outseqs = []
    outints = []
    outdists = []
    outsecs = []
    outnames = []
    outstatists = []
    outcoverage = []
    templates = []
    
    # map rnacompete experiment to unique pdbs
    # accept multiple experiments for one pdb with cutoffs
    # Then choose rnacompete with close to maxmotif (-0.5) and highest identity
    
    
    hasnone = []
    for p, pdb in enumerate(pdbnames):
        loci = np.where(idents[p] >= idcut)[0]
        if len(loci) == 0:
            hasnone.append(pdb)
        elif len(loci) > 0:
            # get the names of the experiments that might work
            rncmpts = np.unique(proteinidname[loci])
            tempmot = tempmotif[p]
            tchain = tempchain[p]
            tempseq = intseq[p]
            tempint = interface[p]
            tempdist = intdistance[p]
            tsempsec = seqstruc[secname == pdb][0]
            tsempseq = secseq[secname == pdb][0]
            maxmotsim = []
            maxcommon = []
            maxchain = []
            maxcover = []
            maxident = []
            alignedseq = []
            alignments = []
            for rncmpt in rncmpts:
                rncmptmot = motifile[motifile[:,0] == rncmpt,1][0]
                sim, alseq, ch = motifsimilarity(rncmptmot, tempmot)
                print alseq
                maxmotsim.append(sim)
                maxcommon.append(alseq)
                # align full rncmpt sequences to pdb and calculate coverage
                
                rncmptseqs = ''.join(fastaseq[proteinidname == rncmpt])
                alignment = pairwise2.align.globalds(rncmptseqs, tempseq, matrix, -22,-2)
                alignments.append(alignment[0])
                rncmptint, iden, cover, rloc= align(alignment[0][0], alignment[0][1], tempint.astype(int), '', 0.4)
                
                rncmptdist, iden, cover, rloc = align(alignment[0][0], alignment[0][1], tempdist,',', 0.4)
                
                alignment = pairwise2.align.globalds(iden,tsempseq, matrix, -22,-2)
                rncmptsec, idenloc, coversec, rlocsec = align(alignment[0][0], alignment[0][1], tsempsec,'', None)
                if len(iden) > 1 and iden != idenloc:
                    print 'Wroing'
                    print iden, len(iden)
                    print idenloc, len(idenloc)
                    print len(rncmptsec)
                    sys.exit()
                    
                maxcover.append(cover)
                alignedseq.append([iden, rncmptint, rncmptdist, rncmptsec])
                maxident.append(100.*float(maxcover[-1])/float(len(iden)))
                maxchain.append([pdb, tchain[ch], rncmpt, str(rloc[0])+'-'+str(rloc[1])])
            # choose experiment with highest coverage and highest motif similarity
            maxmotsim = np.array(maxmotsim)
            maxident = np.array(maxident)
            maxcover = np.array(maxcover)
            finalall = np.where((maxident >= idcut) * (maxmotsim >= motsimcut)* (maxcover >= np.argmax(maxcover)*0.7))[0]
            final = np.where((maxident >= idcut) * (maxmotsim >= max(motsimcut, np.amax(maxmotsim)-0.5)))[0]

            if len(final) == 0:
                hasnone.append(pdb)
            else:
                f = final[np.argmax(maxcover[final])]
                print alignments[f][0]
                print alignments[f][1]
                outobj.write('>'+':'.join(maxchain[f])+':'+','.join(proteinnames[list(proteinnames[:,1]).index(maxchain[f][2]),3:5])+':'+str(maxcover[f])+':'+str(maxident[f])+':'+str(maxmotsim[f])+':'+'-'.join(maxcommon[f])+'\n'+'\n'.join(alignedseq[f])+'\n')
                print '>'+':'.join(maxchain[f])+':'+str(maxcover[f])+':'+str(maxident[f])+':'+str(maxmotsim[f])+':'+'-'.join(maxcommon[f])+'\n'+'\n'.join(alignedseq[f])+'\n'
                print 'Possible', len(finalall)
                for f in finalall:
                    
                    outobj2.write('>'+':'.join(maxchain[f])+':'+','.join(proteinnames[list(proteinnames[:,1]).index(maxchain[f][2]),3:5])+':'+str(maxcover[f])+':'+str(maxident[f])+':'+str(maxmotsim[f])+':'+'-'.join(maxcommon[f])+'\n'+'\n'.join(alignedseq[f])+'\n')
                
    print '\n'.join(hasnone)
    print len(hasnone), len(pdbnames)        
            
            
            
            

       




