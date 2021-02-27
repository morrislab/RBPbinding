import numpy as np
import sys, os
from scipy.spatial import distance
import gzip 
from scipy.stats import pearsonr
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62
from Bio import pairwise2
import operator

pdb, protchains, rnachains = sys.argv[1].split(',')
protchains = [protchains]
rnachains = rnachains.split('-')
scorefile, scorename = sys.argv[2].split(',')
dnum = sys.argv[3]
scorecut = float(sys.argv[4])
if ',' in dnum:
    dnum = np.array(dnum.split(','),dtype = int)
else:
    dnum = [int(dnum)]

dread = False
obj = open(scorefile, 'r').readlines()
for l, line in enumerate(obj):
    if line[0] == '>':
        sname = line[1:].strip()
        if '||' in sname:
            sname = sname.split('||')[1]
        if '__' in sname:
            dr = int(sname[-1])
            sname = sname.split('__')[0]
            dread = True
            
        if sname == scorename:
            if dread:
                if dr in dnum:
                    scseq = obj[l+1].strip()
                    scscore = np.array(obj[l+2].strip().split(','),dtype = float)
                    
            else:
                wseq = obj[l+1].strip()
                if '*' in wseq:
                    loc = np.concatenate([[-1], np.where(np.array(list(wseq)) == '*')[0], [len(wseq)]])
                    sl = loc[dnum[0]]+1
                    el = loc[dnum[-1]+1]
                    scseq = wseq[sl:el]
                else:
                    sl = 0
                    el = len(wseq)
                    scseq = wseq
                scscore = np.array(obj[l+2].strip().split(','),dtype = float)[sl:el]
                if '*' in scseq:
                    wsmask = np.array(list(scseq)) != '*'
                    loc = np.concatenate([[-1], np.where(~wsmask)[0], [len(scseq)]])
                    scseq = ''.join(np.array(list(scseq))[wsmask])
                    for i in range(len(dnum)):
                        # fixes strong edge effects
                        # assigns mimum in center of sequence to edge residues if smaller than minimum
                        lowd = np.amin(scscore[loc[i]+1+5:loc[i+1]-5])
                        lowend = scscore[loc[i]+1:loc[i]+1+6]
                        lowend[lowend < lowd] = lowd
                        scscore[loc[i]+1:loc[i]+1+6] = lowend
                        highend = scscore[loc[i+1]-6:loc[i+1]]
                        highend[highend < lowd] = lowd
                        scscore[loc[i+1]-6:loc[i+1]] = highend
                        if '--normscore' in sys.argv:
                            scscore[loc[i]+1:loc[i+1]] = (scscore[loc[i]+1:loc[i+1]] - np.amin(scscore[loc[i]+1:loc[i+1]]))/np.std(scscore[loc[i]+1:loc[i+1]])
                    scscore = scscore[wsmask]
                else:
                    lowd = np.amin(scscore[6:-5])
                    lowend = scscore[:6]
                    lowend[lowend < lowd] = lowd
                    scscore[:6] = lowend
                    highend = scscore[-5:]
                    highend[highend < lowd] = lowd
                    scscore[-5:] = highend
                    if '--normscore' in sys.argv:
                        scscore = (scscore - np.amin(scscore))/np.std(scscore)
            if '--meanlowscore' in sys.argv:
                scscore = scscore - np.mean(scscore)
                scscore[scscore < 0] = 0.
    

            
def fillgaps(re, ga, ver):
    gax = np.array(ga)-min(ga)
    fullseq=[]
    oo = 0
    for fg in range(max(gax)+1):
        if fg in gax:
            fullseq.append(re[oo])
            oo += 1
        else:
            fullseq.append(ver)
    return np.array(fullseq)

    
# read in pdb sequence
if '.gz' in pdb:
    pdbobj = gzip.open(pdb, 'r')
else:
    pdbobj = open(pdb, 'r')
plines = pdbobj.readlines()

#check for all proteins and RNA chains in the file


modelbarrs = []
header = []
linechain = []
linenum = []
inlines = []
head = True
for i, line in enumerate(plines):
    if line[:4] == 'ATOM' or line[:6] == 'HETATM' or line[:3] == 'TER':
        head = False
        linechain.append(line[21])
        inlines.append(line[:54])
        linenum.append(int(line[6:11]))
    if head:
        header.append(line)
    if line[:5] == 'MODEL':
        modelbarrs.append(i)


        
# if different models exist, only use the first model!!!
if len(modelbarrs) > 1:
    modelbarrs = modelbarrs[:2]
else:
    modelbarrs = [0, len(plines)]
  



dictionary = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'ASN':'N', 'PRO':'P', 'VAL':'V', 'ILE':'I', 'CYS':'C', 'TYR':'Y', 'HIS':'H', 'THR':'T', 'GLY':'G', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'LYS':'K', 'GLN':'Q', 'GLU':'E', 'SER':'S', 'XXX':'X'}

# read in all chains sorted after the chain names
aaseqs = []
aaseqnums = []
aaresids = []
aacoords = []
atomnums = []
atomtypes = []
restarts = []
for k, pchain in enumerate(protchains):
    atomnum = []
    aaseq = []
    aaseqnum = []
    aares = []
    aacoor = []
    aacheck = []
    readout = False
    for p in range(modelbarrs[0], modelbarrs[1]):
        pline = plines[p]
        if pline[:4] == 'ATOM' or pline[:6] == 'HETATM':
            linechain.append(pline[21])
            inlines.append(pline[:54])
            linenum.append(int(pline[6:11]))
            plinecomp = pline.strip().split()
            if pline[21] == pchain:
                readout = True
                if pline[6:11].strip() not in atomnum:
                    Aacid = pline[17:21].strip()
                    if Aacid not in dictionary.keys():
                        #print Aacid, 'Residue not understood'
                        shrot = True
                        skrot = True
                        if pline[17:20] in dictionary.keys():
                            shrot = False
                            Aacid = pline[17:20]
                        if pline[16:19] in dictionary.keys():
                            skrot = False
                            Aacid = pline[16:19]
                        if shrot and skrot:
                            Aacid = 'X'
                    if (pline[12:16].strip() == 'CA' or pline[12:16].strip() == 'C' or pline[12:16].strip() == 'N' or pline[12:16].strip() == 'O') and int(pline[22:26]) not in aaseqnum:
                        aaseq.append(Aacid)
                        aaseqnum.append(int(pline[22:26]))
                    if pline[12:16].strip()[0] != 'H':
                        atomnum.append(pline[6:11].strip())
                        aacoor.append([float(pline[30:38]), float(pline[38:46]), float(pline[46:54])])
                        aares.append(int(pline[22:26]))
                        aacheck.append(pline[12:16].strip())
        if readout and pline[:3] == 'TER':
            linechain.append(pline[21])
            inlines.append(pline[:54])
            linenum.append(int(pline[6:11]))
            break
    restart = np.amin(aaseqnum)
    aaseqnum = np.array(aaseqnum) - restart
    aares = np.array(aares) - restart
    restarts.append(restart)
    atomnums.append(atomnum)
    aaseqs.append(np.array(aaseq))
    aaseqnums.append(aaseqnum)
    aaresids.append(aares)
    aacoords.append(aacoor)
    atomtypes.append(aacheck)
    
# if only main chain is given pdb increase the distance to characterize the interface

if 'CB' not in np.unique(aacheck):
    print '!!! side chains not present in pdb!!! Only coordinates for:', ' '.join(np.unique(aacheck))
    cutoffdist = (np.sqrt(cutoffdist) + 5.)**2
    sys.exit()

nuccoords = []
nucresids = []
nucleotypes = []
nucatomtypes = []
nucatomnum = []
residuesequence = []
residuenumbers = []
nucstarts = []
for k, ligchain in enumerate(rnachains):
    nucoord = []
    natomnum = []
    nucres = []
    nuctypes = []
    nucatoms = []
    residueseq = []
    residuenum = []
    readout = False
    currnuc = -1
    for p, pline in enumerate(plines):
        if pline[:4] == 'ATOM' or pline[:6] == 'HETATM':
            linechain.append(pline[21])
            inlines.append(pline[:54])
            linenum.append(int(pline[6:11]))
            plinecomp = pline.strip().split()
            if pline[21] == ligchain:
                readout = True
                nucoord.append([float(pline[30:38]), float(pline[38:46]), float(pline[46:54])])
                nucres.append(int(pline[22:26]))
                if currnuc != int(pline[22:26]):
                    currnuc = int(pline[22:26])
                    residueseq.append(pline[17:21].strip())
                    residuenum.append(int(pline[22:26]))
                nuctypes.append(pline[17:21].strip())
                natomnum.append(pline[6:11].strip())
                nucatoms.append(pline[12:16].strip())
        if readout and pline[:3] == 'TER':
            linechain.append(pline[21])
            inlines.append(pline[:54])
            linenum.append(int(pline[6:11]))
            break
    nucstart = np.amin(nucres)
    residuenum = np.array(residuenum) - nucstart
    nucstarts.append(nucstart)
    residueseq = np.array(residueseq)
    residuenumbers.append(residuenum)
    residuesequence.append(residueseq)
    nucres = np.array(nucres) - nucstart
    nuccoords.append(nucoord)
    nucresids.append(nucres)
    nucleotypes.append(nuctypes)
    nucatomtypes.append(nucatoms)
    nucatomnum.append(natomnum)

nucatomnum = np.array(nucatomnum)
aacoords = np.array(aacoords)
nuccoords = np.array(nuccoords)
nucleotypes = np.array(nucleotypes)
nucatomtypes = np.array(nucatomtypes)
aaresids = np.array(aaresids)
aaseqs = np.array(aaseqs)
aaseqnums = np.array(aaseqnums)
nucresids = np.array(nucresids)
atomnums = np.array(atomnums)
atomtypes = np.array(atomtypes)

linenum = np.array(linenum)
linechain = np.array(linechain)
inlines = np.array(inlines)




for i, pchain in enumerate(protchains):
    # transform sequence to single character#
    mylen = np.vectorize(len)
    for j in range(len(dictionary)):
        aaseqs[i][aaseqs[i] == dictionary.keys()[j]] = dictionary.values()[j]
    #print len(aaseqs[i]) 
    #print len(scseq), len(scscore)
    if (mylen(aaseqs[i]) > 1).any():
        print 'AA not understood', aaseqs[i]
        sys.exit()
    # align pdb sequence to rncmpt sequence
    alignment = pairwise2.align.localds(''.join(aaseqs[i]), scseq, matrix, -11., -1.)
    pdbscore = -np.ones(len(aaseqs[i]))
    pa = 0
    pb = 0
    sid = 0.
    sclen = 0
    print alignment[0][0]
    print alignment[0][1]
    for t, et in enumerate(alignment[0][0]):
        wt = alignment[0][1][t]
        if wt != '-' and et !='-':
            pdbscore[pa] = scscore[pb]
            sclen +=1
        if et == wt:
            sid += 1.
        if et != '-':
            pa += 1
        if wt != '-':
            pb +=1
    pdbscore[pdbscore != -1] = pdbscore[pdbscore != -1] - np.amin(pdbscore[pdbscore != -1])
    if '--maskpdb' in sys.argv:
        pdbmask = pdbscore != -1
    elif '--maskpdbedges' in sys.argv:
        pdbmask = pdbscore != -1
        pdbedges = np.where(pdbmask)[0][[0,-1]]
        pdbmask[pdbedges[0]:pdbedges[1]] = True
    else:
        pdbmask = np.ones(len(pdbscore)) ==1
    pdbscore[pdbscore == -1] = 0.
    
    pdbbinscore = (pdbscore > scorecut).astype(int)    
    print 'Identity to template', sid/float(len(alignment[0][0])), sid/float(len(aaseqs[i])), sid/float(len(scseq)), sid
    print 'matching scorepositions', sclen
    print 'written', scorename+'_to_'+os.path.split(pdb)[1][:4]+'_'+pchain+'_'+rnachains[0]+os.path.splitext(os.path.split(scorefile)[1])[0]+'.pdb'
    pobj =open(scorename+'_to_'+os.path.split(pdb)[1][:4]+'_'+pchain+'_'+rnachains[0]+os.path.splitext(os.path.split(scorefile)[1])[0]+'.pdb', 'w')
    for hline in header:
        pobj.write(hline)
    curaa = 0
    for ai in range(len(aacoords[i])):
        if aaresids[i][ai] != aaresids[i][max(0, ai-1)]:
            curaa += 1
        if pdbmask[curaa]:
            pobj.write('%54s%6s%6s%12s  \n' % (inlines[linenum == int(atomnums[i][ai])][0], np.around(pdbscore[curaa], 2), np.around(pdbbinscore[curaa],2), atomtypes[i][ai][0]))
    pobj.write('TER %50s%6s%6s%12s  \n' % (inlines[linenum == int(atomnums[i][ai])+1][0][4:54], '', '', ''))
    
    for d, ligchain in enumerate(rnachains):
        curaa = 0
        distmat = distance.cdist(aacoords[i], nuccoords[d], 'euclidean')
        distmask = distmat <= 5.
        rintermap = np.sum(distmask, axis = 0) > 0
        rinterdist = np.around(np.amin(distmat, axis = 0),2)
        for ai in range(len(nuccoords[d])):
            if nucresids[d][ai] != nucresids[d][min(0, ai-1)]:
                curaa += 1
            pobj.write('%54s%6s%6s%12s  \n' % (inlines[linenum == int(nucatomnum[d][ai])][0], rintermap[ai].astype(int), rinterdist[ai], nucatomtypes[d][ai][0]))
        pobj.write('TER %50s%6s%6s%12s  \n' % ('', '', '', '')) #(inlines[linenum == int(nucatomnum[d][ai])+1][0][4:54]
            

