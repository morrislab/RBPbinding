import numpy as np
import sys, os
from scipy.spatial import distance
import gzip 

pdb = sys.argv[1]

if '--interdistcut' in sys.argv:
    cutint = True
    cutoffdist = float(sys.argv[sys.argv.index('--interdistcut')+1])**2
    outcut = '_lt'+sys.argv[sys.argv.index('--interdistcut')+1]
elif '--interdist' in sys.argv:
    cutint = False
    outcut = '_interface_dist'
else:
    cutint = True
    cutoffdist = 25.
    outcut = '_lt5'
    
    
# only use that option if you have one RNA and one protein interacting !!!!!!
pernuc = False
if '--pernucleotide' in sys.argv:
  pernuc = True

prchain = None
if '--proteinchains' in sys.argv:
    prchain = sys.argv[sys.argv.index('--proteinchains')+1]

outdir = ''
if '--outdir' in sys.argv:
    outdir = sys.argv[sys.argv.index('--outdir')+1]
    if outdir[-1] != '/':
        outdir = outdir+'/'

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
protchains = []
rnachains = []
allchains = []
controlchains = []
controlseq = []
modelbarrs = []

header = []
linechain = []
linenum = []
inlines = []
head = True
for i, line in enumerate(plines):
    if line[:4] == 'ATOM' or line[:6] == 'HETATM' or line[:3] == 'TER':
        head = False
    if head:
        header.append(line)

    
    if line.split()[0] == 'DBREF':
        line = line.split()
        allchains.append([line[2], line[6]+'_'+line[3]+'-'+line[4]])
    if line[:6] == 'SEQRES':
        controlchains.append(line.split()[2])
        controlseq.append(line.strip().split()[4:])
    if line[:5] == 'MODEL':
        modelbarrs.append(i)

controlseq = np.array(controlseq)
controlchains = np.array(controlchains)


for p, potchain in enumerate(np.unique(controlchains)):
    cseqs = np.concatenate(controlseq[controlchains == potchain[0]])
    if np.isin(cseqs, ['A', 'C', 'G', 'U']).any():
        rnachains.append([potchain, pdb[:4]])
    else:
        protchains.append([potchain, pdb[:4]])
if len(protchains) == 0:
    print 'No protein chains found' 
    sys.exit()
if len(rnachains) == 0:
    print 'No RNA chain found'
    sys.exit()

if '--rnachains' in sys.argv:
    rchain = sys.argv[sys.argv.index('--rnachains')+1]
    rnaconst = sys.argv[sys.argv.index('--rnachains')+2].split(',')
    crnachain = '_R'+rchain+'nt'+rnaconst[0]+'-'+rnaconst[1]
else:
    rchain = None
    crnachain = ''
    
# if different models exist, only use the first model!!!
if len(modelbarrs) > 1:
    modelbarrs = modelbarrs[:2]
else:
    modelbarrs = [0, len(plines)]
  
    

protchains = np.array(protchains)

unprotchains = np.unique(protchains[:,0])
if len(unprotchains) < len(protchains):
    newprotchains = []
    for unprotchain in unprotchains:
        protname = ''
        for protchain in protchains:
            if protchain[0] == unprotchain:
                protname += protchain[1]+'_' 
        newprotchains.append([unprotchain, protname[:-1]])
    protchains = np.array(newprotchains)




dictionary = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'ASN':'N', 'PRO':'P', 'VAL':'V', 'ILE':'I', 'CYS':'C', 'TYR':'Y', 'HIS':'H', 'THR':'T', 'GLY':'G', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'LYS':'K', 'GLN':'Q', 'GLU':'E', 'SER':'S', 'XXX':'X'}


print protchains, rnachains

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
            if pline[21] == pchain[0]:
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
    currnuc = -100
    for p, pline in enumerate(plines):
        if pline[:4] == 'ATOM' or pline[:6] == 'HETATM':
            linechain.append(pline[21])
            inlines.append(pline[:54])
            linenum.append(int(pline[6:11]))
            plinecomp = pline.strip().split()
            if pline[21] == ligchain[0]:
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
    residuenum = np.array(residuenum) # - nucstart
    nucres = np.array(nucres)
    nucstarts.append(nucstart)
    nnucres = -np.ones(len(nucres))*1000
    for n, nres in enumerate(residuenum):
        nnucres[nucres == nres] = n
    residueseq = np.array(residueseq)
    residuenumbers.append(np.arange(len(residuenum), dtype = int)) #residuenum)
    residuesequence.append(residueseq)
    nucres = np.array(nnucres, dtype = int) #np.array(nucres) - nucstart
    nuccoords.append(np.array(nucoord))
    nucresids.append(np.array(nucres))
    nucleotypes.append(np.array(nuctypes))
    nucatomtypes.append(np.array(nucatoms))
    nucatomnum.append(np.array(natomnum))


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
    if prchain is None or prchain == pchain[0]:
        mylen = np.vectorize(len)
        for j in range(len(dictionary)):
            aaseqs[i][aaseqs[i] == dictionary.keys()[j]] = dictionary.values()[j]
            
        if (mylen(aaseqs[i]) > 1).any():
            print 'AA not understood', aaseqs[i]
            sys.exit()
        hascontact = [False]
        distmats = []
        for ii in range(len(rnachains)):    
            # calculate distance between all protein and rna molecules
            if rchain is not None:
                takerna = False
                if rchain.isdigit():
                    if ii == int(rchain):
                        takerna = True
                else:
                    if rchain == rnachains[ii][0]:
                        takerna = True
                if takerna:
                    constresidues = np.arange(int(rnaconst[0]), int(rnaconst[1]), dtype = int)
                    nucresmask = np.isin(nucresids[ii], constresidues)
                    nuccoordsd = nuccoords[ii][nucresmask]
                else:
                    nuccoordsd = nuccoords[ii]
            else:
                takerna = True
                nuccoordsd = nuccoords[ii]
            distmat = distance.cdist(aacoords[i], nuccoordsd, 'sqeuclidean')
            ### check if RNA and protein have atoms in contact
            if np.sum(distmat <= cutoffdist) > 1 and takerna:
                hascontact.append(True)
            distmats.append(distmat)
        keeprna = []
        if np.array(hascontact).any():
            wobj =open(outdir+pchain[1]+'_'+pchain[0]+crnachain+outcut+'.int', 'w')
            fobj =open(outdir+pchain[1]+'_'+pchain[0]+crnachain+'.fasta', 'w')
            pobj =open(outdir+pchain[1]+'_'+pchain[0]+crnachain+outcut+'-interface.pdb', 'w')
            for hline in header:
                #print hline
                pobj.write(hline)
            interfacemap = np.zeros(np.shape(aaseqs[i]), dtype = int)
            interfacedistmap = np.ones(np.shape(aaseqs[i]), dtype = float)*100.
            interfacebackmap = np.chararray(np.shape(aaseqs[i]), itemsize = 1)
            interfacebackmap[:] = '*'
            interfaceresmap = np.chararray(np.shape(aaseqs[i]), itemsize = 1)
            interfaceresmap[:] = '*'
            interfaceresrnamap = np.chararray(np.shape(aaseqs[i]), itemsize = 28)
            interfaceresrnamap[:] = '*'
            rnasequences = []
            if pernuc == False:
                newresiduenumbers = []
                for d, distmat in enumerate(distmats):
        
                    if rchain is not None:
                        takerna = False
                        if rchain.isdigit():
                            if d == int(rchain):
                                takerna = True
                        else:
                            if rchain == rnachains[d][0]:
                                takerna = True
                        if takerna:
                            constresidues = np.arange(int(rnaconst[0]), int(rnaconst[1]), dtype = int)
                            residuenumbermask = np.isin(residuenumbers[d], constresidues)
                            residuenumbersd = residuenumbers[d][residuenumbermask]
                            residuesequenced = residuesequence[d][residuenumbermask]
                            nucresmask = np.isin(nucresids[d], constresidues)
                            nucresidsd = nucresids[d][nucresmask]
                            nuccoordsd = nuccoords[d][nucresmask]
                            nucleotypesd = nucleotypes[d][nucresmask]
                            nucatomtypesd = nucatomtypes[d][nucresmask]
                            nucatomnumd = nucatomnum[d][nucresmask]
                            distmat = distance.cdist(aacoords[i], nuccoordsd, 'sqeuclidean')
                    else:
                        takerna = True
                        residuenumbersd = residuenumbers[d]
                        residuesequenced = residuesequence[d]
                        nucresidsd = nucresids[d]
                        nuccoordsd = nuccoords[d]
                        nucleotypesd = nucleotypes[d]
                        nucatomtypesd = nucatomtypes[d]
                        nucatomnumd = nucatomnum[d]
                        
                    distmask = distmat <= cutoffdist
                    if np.sum(distmask) > 1 and takerna:
                                
                        rintermap = np.sum(distmask, axis = 0) > 0
                        #print len(rnaseq)
                        fullrnaseq =  '-'.join(residuesequenced)
                        
                    
                        
                        
                        rnasugar = np.isin(nucatomtypesd, ["O5'", "C5'","O4'","C4'","O3'","C3'","O2'","C2'"])
                        rnaphosphor = np.isin(nucatomtypesd, ["OP1", "P", "OP2"])
                        # interactions with sugar
                        suginteraction = np.sum(distmask[:, rnasugar], axis = 1)
                        # interactions with phosphor
                        phosinteraction = np.sum(distmask[:, rnaphosphor], axis = 1)
                        allinteraction = np.sum(distmask, axis = 1)
                        # base-interactions/ base interactions shouldn't be excluded if sugar or phosphor exists
                        diffinteraction = allinteraction - suginteraction - phosinteraction
                        suginteraction = suginteraction > 0
                        phosinteraction = phosinteraction > 0 
                        diffinteraction = diffinteraction > 0 
                        suginteraction = np.unique(np.array(aaresids[i])[suginteraction])
                        suginteraction = np.isin(aaseqnums[i], suginteraction)
                        phosinteraction = np.unique(np.array(aaresids[i])[phosinteraction])
                        phosinteraction = np.isin(aaseqnums[i], phosinteraction)
                        
                        rnaseq = np.chararray(len(residuenumbersd))
                        resnumstart = np.amin(residuenumbersd)
                        residueindexes = list(residuenumbersd)
                        rnaseq[:] = 'N'
                        rnabase = ~np.isin(nucatomtypesd, ["O5'", "C5'","O4'","C4'","O3'","C3'","O2'","C2'","OP1", "P", "OP2"])
                        basenames = []
                        for r, residn in enumerate(nucresidsd):
                            basenames.append(rnachains[d][0]+'.'+nucleotypesd[r]+'^'+str(residn))
                        basenames = np.array(basenames)
                        rnaassign = np.chararray(len(aaseqnums[i]), itemsize = 28)
                        rnaassign[:] = '*'
                        for a, aaseqn in enumerate(aaseqnums[i]):
                            aaseqatm = aaresids[i] == aaseqn
                            aaresatm = np.sum(distmask[aaseqatm], axis = 0) > 0
                            if np.sum(aaresatm) > 0:
                                leftbases = aaresatm*rnabase
                                if np.sum(leftbases) > 0:
                                    boundrnabase = np.unique(basenames[leftbases])
                                    for basen in boundrnabase:
                                        rnaseq[residuenumbersd == int(basen.split('^')[-1])] = basen.split('^')[0].split('.')[1]
                                    boundrnabase = '-'.join(boundrnabase)
                                    rnaassign[a] = boundrnabase
                                    
                        # make shorter/ easier to visualize version of rnaseq
                        motiflocations = []
                        strna = None
                        endrna = None
                        countN = 0
                        countseq = 0
                        startfound = True
                        countNend = 0
                        countseqend = 0
                        endfound = True
                        posspec = []
                        for n in range(len(rnaseq)):
                            if strna is None and rnaseq[n] != 'N':
                                strna = n
                            if endrna is None and rnaseq[-1-n] != 'N':
                                endrna = len(rnaseq)-n
                        rnaseq = rnaseq[strna:endrna]
                        if strna is None:
                            motifstarts = 0
                        else:
                            motifstarts = strna
                        if endrna is None:
                            motifends = len(rnaseq)
                        else:
                            motifends = endrna
                        strna = 0
                        for n in range(len(rnaseq)):
                            if rnaseq[n] == 'N':
                                countN += 1
                            if countN > 0 and rnaseq[n] != 'N':
                                if countN > 4:
                                    posspec.append(''.join(np.array(rnaseq[strna:strna+countseq-countN])))
                                    motiflocations.append(str(resnumstart+motifstarts+strna)+'-'+str(resnumstart+motifstarts+strna+countseq-countN))
                                    strna = n
                                    countseq = 0
                                countN = 0
                            countseq += 1
                            if n == len(rnaseq)-1 and countseq > 0:
                                posspec.append(''.join(np.array(rnaseq[strna:strna+countseq])))
                                motiflocations.append(str(resnumstart+motifstarts+strna)+'-'+str(resnumstart+motifstarts+strna+countseq))
                        include = False
                        for pmot in posspec:
                            if len(pmot)>=3:
                                include = True
                        if include:
                            keeprna.append(d)
                            rnasequences.append('### '+str(rnachains[d][0])+' '+fullrnaseq)
                            rnasequences.append('## '+','.join(posspec))
                            rnasequences.append('# '+','.join(motiflocations))
                            
                            diffinteraction = np.unique(np.array(aaresids[i])[diffinteraction])
                            diffinteraction = np.isin(aaseqnums[i], diffinteraction)
                            
                            intermap = np.sum(distmask, axis = 1) > 0
                            inbackbone = np.isin(atomtypes[i], ['CA', 'N', 'O', 'C'])
                            
                            # protein backbone interactions
                            interbackbone = inbackbone & intermap
                            interresidue = np.invert(inbackbone) & intermap

                            interfaceresidues = np.unique(np.array(aaresids[i])[intermap])
                            intermask =  np.isin(aaseqnums[i], interfaceresidues) 
                            
                            interfacebresidues = np.unique(np.array(aaresids[i])[interbackbone > 0.])
                            interBmask =  np.isin(aaseqnums[i], interfacebresidues)
                            
                            interfaceRresidues = np.unique(np.array(aaresids[i])[interresidue])
                            interRmask =  np.isin(aaseqnums[i], interfaceRresidues)
                            
                            interfacemap[intermask] = 1
                            interfacebackmap[interBmask] = 'B'
                            interfaceresmap[interRmask] = 'R'
                            interfaceresrnamap[suginteraction] = 'S'
                            interfaceresrnamap[phosinteraction] = 'P'
                            
                            interfaceresrnamap[diffinteraction] = rnaassign[diffinteraction] # 'N'
                            
                            interdist = np.amin(distmat, axis = 1)
                            for ai in aaseqnums[i]:
                                ###NEEDED fix
                                interfacedistmap[aaseqnums[i] == ai] = min(interfacedistmap[aaseqnums[i] == ai], np.around(np.sqrt(np.amin(interdist[aaresids[i] == ai])),2))
                
                if len(rnasequences) > 0:
                    wobj.write('\n'.join(np.array(rnasequences))+'\n')
                    wobj.write('>'+pchain[1]+'_'+pchain[0]+'\n')
                    wobj.write(','.join(aaseqs[i])+'\n')
                    wobj.write(','.join(interfacemap.astype(str))+'\n')
                    wobj.write(','.join(interfaceresmap.astype(str))+'\n')
                    wobj.write(','.join(interfacebackmap.astype(str))+'\n')
                    wobj.write(','.join(interfaceresrnamap.astype(str))+'\n')
                    wobj.write(','.join(interfacedistmap.astype(str))+'\n')
                    fobj.write('>'+pchain[1]+'_'+pchain[0]+'\n'+''.join(aaseqs[i])+'\n')
            
            # generate pdb with the interface distances and the rna sequences that cover enough residues
            print pchain, keeprna, rnasequences 
            curaa = 0
            for ai in range(len(aacoords[i])):
                if aaresids[i][ai] != aaresids[i][max(0, ai-1)]:
                    curaa += 1
                pobj.write('%54s%6s%6s%12s  \n' % (inlines[linenum == int(atomnums[i][ai])][0], np.around(interfacemap[curaa], 2), np.around(interfacedistmap[curaa],2), atomtypes[i][ai][0]))
            pobj.write('TER %50s%6s%6s%12s  \n' % (inlines[linenum == int(atomnums[i][ai])+1][0][4:54], '', '', ''))
            #pobj.write('TER %7s  %-4s%3s%2s%4s    %8s%8s%8s%6s%6s%12s  \n' % (int(atomnums[i][ai])+1, '', dictionary.keys()[dictionary.values().index(aaseqs[i][curaa])], pchain[0], aaseqnums[i][curaa]+restarts[i], '','','','','',''))    
            nrn = -1
            
            for d in keeprna:
                distmat = distmats[d]
                
                if rchain is not None:
                    takerna = False
                    if rchain.isdigit():
                        if d == int(rchain):
                            takerna = True
                    else:
                        if rchain == rnachains[d][0]:
                            takerna = True
                    if takerna:
                        constresidues = np.arange(int(rnaconst[0]), int(rnaconst[1]), dtype = int)
                        nucresmask = np.isin(nucresids[d], constresidues)
                        nucresidsd = nucresids[d][nucresmask]
                        nuccoordsd = nuccoords[d][nucresmask]
                        nucatomtypesd = nucatomtypes[d][nucresmask]
                        nucatomnumd = nucatomnum[d][nucresmask]
                        distmat = distance.cdist(aacoords[i], nuccoordsd, 'sqeuclidean')
                else:
                    takerna = True
                    nucresidsd = nucresids[d]
                    nuccoordsd = nuccoords[d]
                    nucatomtypesd = nucatomtypes[d]
                    nucatomnumd = nucatomnum[d]
            
                
                distmask = distmat <= cutoffdist
                if np.sum(distmask) > 1 and takerna:
                    nrn += 1
                    rintermap = np.sum(distmask, axis = 0) > 0
                    rinterdist = np.around(np.sqrt(np.amin(distmat, axis = 0)),2)
                    curaa = 0
                    #print len(nuccoords[d]), len(newresiduenumbers[nrn]), len(nucresids[d]), len(inlines), len(rinterdist), len(rintermap), len(nucatomtypes[d]), len(atomnums[i])
                    for ai in range(len(nuccoordsd)):
                        #if nucresids[d][ai] != nucresids[d][min(0, ai-1)]:
                        if nucresidsd[ai] != nucresidsd[min(0, ai-1)]:
                            curaa += 1
                        #print len(nucatomnum[d])#[ai])
                        pobj.write('%22s%4s%28s%6s%6s%12s  \n' % (inlines[linenum == int(nucatomnumd[ai])][0][:22], nucresidsd[ai], inlines[linenum == int(nucatomnumd[ai])][0][26:54], rintermap[ai].astype(int), rinterdist[ai], nucatomtypesd[ai][0]))
                    pobj.write('TER %18s%4s%28s%6s%6s%12s  \n' % (inlines[linenum == int(nucatomnumd[ai])][0][4:22], nucresidsd[ai], inlines[linenum == int(nucatomnumd[ai])][0][26:54], '', '', ''))
                    

            wobj.close()
            fobj.close()
