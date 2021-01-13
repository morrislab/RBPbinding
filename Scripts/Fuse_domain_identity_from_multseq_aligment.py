#fuse_domain_identity.py
import numpy as np
import sys, os


if '--multiprocessing' in sys.argv:
    from joblib import Parallel, delayed
    import multiprocessing    
    mpross = True
    if '--numcores' in sys.argv:
        num_cores = int(sys.argv[sys.argv.index('--numcores')+1])
    else:
        num_cores = multiprocessing.cpu_count() 
else:
    mpross = False
    

   


# creates a matrix based on average domain identity, parwise alignment, domain-ident and structural construction

if '--domainidentities' in sys.argv:
    domainidents = sys.argv[sys.argv.index('--domainidentities') +1]
    dfiletype = os.path.splitext(domainidents)[1]
else:
    print 'Domain identities must be given'
    sys.exit()
    
#@jit(nopython = True, parallel = True)
def add_pair(arr, pair):
    outarr = []
    for ar in arr:
        for pi in pair:
            apend = True
            for arsub in ar:
                if pi[0] == arsub[0] or pi[1] == arsub[1]:
                    apend = False
            if apend == True:
                outarr.append(np.append(ar,[pi], axis =0 ))
    return outarr

#@jit(nopython = True, parallel = True)
def highestfunc(mat, n, lenweight):
        ouid = []
        pairs = []
        for i in range(len(mat)):
            for j in range(len(mat[0])):
                pairs.append([i,j])
        
        for pair in pairs:
            ouid.append([pair])
        
        for k in range(min(min(np.shape(mat)),n)-1):
            ouid = add_pair(ouid, pairs)

        outval = [0.]
        for oui in ouid:

            idou = 0.
            for pou in oui:
                idou += mat[pou[0],pou[1]]
            outval.append(idou)
        shortside = min(len(mat), len(mat[0]))    
        n = min(n, shortside)
        return max(outval)/n

#@jit(nopython = True, parallel = True)
def highfunc(mat, n, lenweight):
        ouid = []
        shpe = (len(mat), len(mat[0]))
        shortlen = min(shpe)
        longlen = max(shpe)
        sideord = np.argmin(shpe)
        if sideord == 1:
            mat = mat.T
            lenweight = lenweight.T
        for o in range(longlen-shortlen+1):
            identchain = mat[0,o]*lenweight[0,o]
            inorm = np.copy(lenweight[0,o])
            for p in range(1, shortlen):
                identchain += mat[p,o+p]*lenweight[p,o+p]
                inorm += lenweight[p,o+p]
            ouid.append(identchain/max(n,inorm))  
        return max(ouid)


#### read in files with proteins and its domains and the domain-domain similarities

#print plist.domlen
#sys.exit()

if dfiletype == '.npz':
    loadmat = np.load(domainidents)
    didents = loadmat['identmat']
    dnames = loadmat['names']
    filenames = loadmat.files
    if 'names2' in filenames:
        dnames2 = loadmat['names2']
        dnames2 = np.array(dnames2)
        #print dnames, dnames2
        if len(np.intersect1d(dnames, dnames2)) == 0:
            print 'Axis have no domain in common'
            two_set = True
            didents = didents.T
            print np.shape(didents)
            print len(dnames), len(dnames2)
        else:
            two_set = False

    loadmat = []
else:
    two_set = False
    dnames = open(domainidents, 'r').readline().strip().split()[1:]
    didents = np.genfromtxt(domainidents)

dnames = np.array(dnames)

# See if every domain only occurs once
uniquedomains = np.unique(dnames, return_index = True, return_counts = True)

if len(uniquedomains[0]) != len(dnames):
    print 'Attention!!!! Some domains occurred twice!!!'
    dnames = uniquedomains[0]
    print uniquedomains
    print dnames[uniquedomains[2] > 1]
    didents = didents[uniquedomains[1]] 
    if two_set == False:
        didents = didents[:, uniquedomains[1]]

if '--proteinlist' in sys.argv:
    protlistfile = sys.argv[sys.argv.index('--proteinlist')+1]
    psubset = int(sys.argv[sys.argv.index('--proteinlist')+2])
    setname = '_set'+str(psubset)
    protlistobj = open(protlistfile, 'r').readlines()
    plist = protlistobj[0].strip().split()
    psublist = protlistobj[psubset].strip().split()
    plistoffset = list(plist).index(psublist[0])
    proteinnames = []
    if '__' in dnames[0]:
        for d, dname in enumerate(dnames):
            proteinnames.append(dname.rsplit('__',1)[0])
    elif dnames[0].rsplit('_',2)[0] in plist:
        for d, dname in enumerate(dnames):
            dname = dname.rsplit('_',2)
            proteinnames.append(dname[0])
    elif dnames[0].rsplit('_',1)[0] in plist:
        for d, dname in enumerate(dnames):
            dname = dname.rsplit('_',1)
            proteinnames.append(dname[0])
    else:
        print dnames[0], plist[0], 'not compatible'
        sys.exit()

elif two_set:
    setname='2sets'
    proteinnames = []
    if '__' in dnames[0]:
        for d, dname in enumerate(dnames):
            proteinnames.append(dname.rsplit('__',1)[0])
    else:
        for d, dname in enumerate(dnames):
            dname = dname.rsplit('_',2)
            proteinnames.append(dname[0])
    plist = np.unique(proteinnames)
    
    proteinnames2 = []
    if '__' in dnames2[0]:
        for d, dname in enumerate(dnames2):
            proteinnames2.append(dname.rsplit('__',1)[0])
    else:
        for d, dname in enumerate(dnames2):
            dname = dname.rsplit('_',2)
            proteinnames2.append(dname2[0])
    plist = np.unique(proteinnames)
    psublist = np.unique(proteinnames2)
    plistoffset = 0
    
else:
    setname = ''
    proteinnames = []
    if '__' in dnames[0]:
        for d, dname in enumerate(dnames):
            proteinnames.append(dname.rsplit('__',1)[0])
    else:
        for d, dname in enumerate(dnames):
            dname = dname.rsplit('_',2)
            proteinnames.append(dname[0])
    plist = np.unique(proteinnames)
    psublist = np.copy(plist)
    plistoffset = 0

         

proteinnames = np.array(proteinnames)
domainnumbers = []
pldomains = []
for pl, pli in enumerate(plist):
    domainnumbers.append(len(np.where(proteinnames == pli)[0]))
    domainlocation = np.where(proteinnames == pli)[0]
    pldomains.append(domainlocation[np.argsort(dnames[domainlocation])])
domainnumbers = np.array(domainnumbers)
totprots = len(plist)
pldomains = np.array(pldomains)

if two_set:
    proteinnames2 = np.array(proteinnames2)
    domainnumbers2 = []
    pldomains2 = []
    for pl, pli in enumerate(psublist):
        domainnumbers2.append(len(np.where(proteinnames2 == pli)[0]))
        domainlocation = np.where(proteinnames2 == pli)[0]
        pldomains2.append(domainlocation[np.argsort(dnames2[domainlocation])])
    domainnumbers2 = np.array(domainnumbers2)
    totprots2 = len(psublist)
    pldomains2 = np.array(pldomains2) 

else:
    domainnumbers2 = np.copy(domainnumbers)
    totprots2 = totprots
    pldomains2 = np.copy(pldomains)

#print pldomains2, domainnumbers2, totprots2, pldomains, domainnumbers, totprots
#sys.exit()

if '--lendoms' in sys.argv:
    infofile = sys.argv[sys.argv.index('--lendoms')+1]
    iobj = open(infofile, 'r')
    ipnames = []
    ipdoms = []
    domlengths = []
    for doolen in domainnumbers:
        domlengths.append(np.ones(doolen))
    np.array(domlengths)
    for i, iline in enumerate(iobj):
        iline = iline.strip().split()
        doolen = []
        init = False
        if iline[0] in plist:
            ipname = iline[0]
            init = True
        elif iline[0][:-1] in plist:
            init = True
            ipname = iline[0][:-1]
        elif iline[0][1:-1] in plist:
            init = True
            ipname = iline[0][1:-1]
        else:
            print 'Cannot find', iline[0], 'in idfile'
        if init:
            pid = list(plist).index(ipname)
            pdlen = domainnumbers[pid]
            for j in range(pdlen):
                domlengths[pid][j] = len(iline[-pdlen+j])
        
else:
    domlengths = []
    for doolen in domainnumbers:
        domlengths.append(np.ones(doolen))
    
    if two_set:
        domlengths2 = []
        for doolen in domainnumbers2:
            domlengths2.append(np.ones(doolen))
        domlengths2 = np.array(domlengths2)
    else:
        domlengths2 = np.array(domlengths)

domlengths = np.array(domlengths)        



#dnames = list(dnames)
multfac = 100.
if np.amax(didents[0]) > 1.01:
    multfac = 1.
didents = didents * multfac




###### start main program

if '--ident_option' in sys.argv:
    # highest domain identity, two highest domains, two highest in order, all domain average without order, all domain in order.
    identopt = sys.argv[sys.argv.index('--ident_option')+1]
    outprot = np.zeros((totprots, len(psublist)))
else:
	print "Define what you mean by domain identity"
	sys.exit()

if "--outname" in sys.argv:
    outname = sys.argv[sys.argv.index("--outname")+1]+setname
else:
    outname =os.path.splitext(domainidents)[0]+'_combined_'+identopt+setname+'.txt'



if identopt == "highdom":
    #@jit(nopython = True)#, parallel = True)
    def maxdomid(domains2, iy):
        if two_set:
            startsec = 0
        else:
            startsec = plistoffset+iy+1
        hidentarray = []
        for j in range(startsec, totprots):
            domains1 = pldomains[j]
            highident = []
            for k, dom1 in enumerate(domains1):
                for l, dom2 in enumerate(domains2):
                    cident =  didents[dom1,dom2]
                    highident.append(cident)
            hidentarray.append(np.amax(highident))
        return hidentarray, iy, startsec
    
    if mpross:
        results = Parallel(n_jobs=num_cores)(delayed(maxdomid)(pldomains2[plistoffset+i], i) for i in range(len(psublist)))
    else:
        results = []
        for i in range(len(psublist)):
            results.append(maxdomid(pldomains2[plistoffset+i], i))
    
    #print len(results), len(results[0]), len(results[1])
    for result, ri, si in results:
        outprot[si:, ri] = result
        


#order compares the domains in order 1,2 to 1,2 or 2,3 for example
elif identopt == "twohigh":
    orden = False
    def numdomain(do1, do2):
        return 2
elif identopt == "twohighorder":
    orden = True
    def numdomain(do1, do2):
        return 2
elif identopt == "alldomain":
    orden = False
    def numdomain(do1, do2):
        return max(do1, do2)
elif identopt == "alldomainorder":
    orden = True
    def numdomain(do1, do2):
        return max(do1, do2)
elif identopt == "shortprot":
    orden = False
    def numdomain(do1, do2):
        return min(do1, do2)
elif identopt == "shortprotorder":
    orden = True
    def numdomain(do1, do2):
        return min(do1, do2)
else:
    print identopt, 'not understood'
    sys.exit()

if identopt != "highdom":
    def mixidentfunc(domains2, dlen2, iy):
        if two_set:
            startsec = 0
        else:
            startsec = plistoffset+iy+1
        hidentarray = []
        for j in range(startsec, totprots):
            domains1 = pldomains[j]
            dlen1 = domlengths[j]
            highident = 0.
            idmat = np.zeros((len(domains1), len(domains2)))
            lenmat = np.zeros((len(domains1), len(domains2)))
            for k, dom1 in enumerate(domains1):
                for l, dom2 in enumerate(domains2):
                    #print dom1, dom2
                    iden = didents[dom1,dom2]
                    idmat[k,l] = iden
                    lenmat[k,l] = max(dlen1[k], dlen2[l])
            ne = numdomain(np.sum(dlen1), np.sum(dlen2))
            if orden:
                highident = highfunc(idmat, ne, lenmat)
            else:
                highident = highestfunc(idmat, ne, lenmat)
            #print highident, plist[j]
            hidentarray.append(highident)
        #sys.exit()
        return hidentarray, iy, startsec
    
    if mpross:
        results = Parallel(n_jobs=num_cores)(delayed(mixidentfunc)(pldomains2[plistoffset+i], domlengths2[plistoffset+i], i) for i in range(len(psublist)) )
    else:
        results = []
        for i in range(len(psublist)):
            #print psublist[i]
            results.append(mixidentfunc(pldomains2[plistoffset+i], domlengths2[plistoffset+i], i))
            #print plist[np.argmax(results[-1][0])], np.amax(results[-1][0])
            
    for result, ri, si in results:
        #print np.shape(outprot)
        #print len(plist), len(psublist), len(result), ri, si
        outprot[si:, ri] = result
    

    
shapeout = np.shape(outprot)
if shapeout[0] == shapeout[1]:
    for i in range(len(psublist)):
        for j in range(plistoffset+i+1, totprots):
            outprot[i,j] = outprot[j,i]
                
np.fill_diagonal(outprot, 100.)

if '--savenpz' in sys.argv:
    print os.path.splitext(outname)[0]+'.npz'
    np.savez_compressed(os.path.splitext(outname)[0]+'.npz', identmat = outprot, names = plist, subnames = psublist)
else:
    print 'Saved as', os.path.splitext(outname)[0]+'.txt'
    np.savetxt(os.path.splitext(outname)[0]+'.txt', outprot, fmt = "%3.2f", header = " ".join(psublist), footer = " ".join(plist) )

if '--savemaxpair' in sys.argv:
    maxobj = open(os.path.splitext(outname)[0]+'_maxstats.txt', 'w')
    for p in range(len(psublist)):
        maxobj.write(psublist[p].split('|',1)[0]+' '+ plist[np.argmax(outprot[:,p])]+' '+str(np.amax(outprot[:,p]))+'\n')





