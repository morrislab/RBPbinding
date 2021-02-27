#cisbpstats.py
import numpy as np
import sys, os
import glob

species = np.sort(glob.glob('*/'))

pmeasobj = np.genfromtxt('../S1.In_Progress.v5.For_all.DR-2_names.txt', dtype = str)
    

dobj = np.genfromtxt('Domaintypes.txt', dtype = str)
domclass = []
for d, dline in enumerate(dobj):
    domclass.append([dline[-1].split('-')[0]])
dobj = np.append(dobj,domclass, axis = 1)

# Similarity file
simfile = sys.argv[1]

# highestorigfile
highfile = sys.argv[2]

# Identity file
idfile = sys.argv[3]

# Name file
namefile = sys.argv[4]
# similarity cutoff: 0.127
simcut = float(sys.argv[5])
# Identity cutoff
# 70 and 50
idcut = float(sys.argv[6])
# Motif clustering
clustfile = sys.argv[7]

outname = sys.argv[8]

mclust = np.load(clustfile)
clustnames = mclust['names'].astype(str)
clusters = mclust['clusters']

if '--maxmeasured' in sys.argv:
    maxexp = int(sys.argv[sys.argv.index('--maxmeasured')+1])
    keep = []
    for p, pmeas in enumerate(pmeasobj):
        if int(pmeas[1][-4:]) <= maxexp:
            keep.append(p)
    pmeasobj = pmeasobj[keep]



outname = 'CisBP_reconstruction-stats-jple'+str(simcut)+'-id'+str(idcut)+'-'+outname+'-'+os.path.splitext(os.path.split(clustfile)[1])[0]+'_nspecies'+str(len(species))+'.csv'

wobj = open(outname, 'w')

wobj.write('Species, Nrbpprot, NprotRRM, NprotKH, Nmeas, NmeasRRM, NmeasKH, Nrectot, Nrecjple, Nrecid, Nrecinter, NrectotRRM, NrecjpleRRM, NrecidRRM, NrecinterRRM, NrectotKH, NrecjpleKH, NrecidKH, NrecinterKH, Nrncmptjple, Nrncmptid, Nmotifjple, Nmotifid\n')

for s, specy in enumerate(species):
    if os.path.isfile(specy+simfile) and os.path.isfile(specy+specy[:-1]+idfile) and os.path.isfile(specy+specy[:-1]+namefile):
        print specy
        sfile = np.genfromtxt(specy+simfile, dtype = str)
        hfile = np.genfromtxt(specy+highfile, dtype = str)
        ifile = np.genfromtxt(specy+specy[:-1]+idfile, dtype = str)
        pfile = np.genfromtxt(specy+specy[:-1]+namefile, dtype = str, delimiter = '', comments = None)
        
        
        rbpprot = np.unique(pfile[:,1])
        Nrbpprot = len(rbpprot)
        
        khprots = np.unique(dobj[np.isin(dobj[:,1], rbpprot) * (dobj[:,-1] == 'KH'),1])
        khprotsp = dobj[np.isin(dobj[:,1], rbpprot) * (dobj[:,-1] == 'KH'),0]
        NprotKH = len(khprots)
        
        rrmprots = np.unique(dobj[np.isin(dobj[:,1], rbpprot) * (dobj[:,-1] == 'RRM'),1])
        rrmprotsp = dobj[np.isin(dobj[:,1], rbpprot) * (dobj[:,-1] == 'RRM'),0]
        NprotRRM = len(rrmprots)
        
        measured, mesindex = np.unique(pmeasobj[pmeasobj[:, -2] == specy.strip('/'),0], return_index = True)
        Nmeas = len(measured)
        
        NmeasKH = int(np.sum(np.isin(measured, khprots)))
        NmeasRRM = int(np.sum(np.isin(measured, rrmprots)))
        
        recid = ifile[ifile[:, -1].astype(float) >= idcut,0]
        recjmf = sfile[sfile[:, -1].astype(float) <= simcut,0] 
        recidjmf = np.union1d(recid, recjmf)
        recidjmfint = np.intersect1d(recid, recjmf)
        
        Nrectot = len(recidjmf)
        Nrecjmf = len(recjmf)
        Nrecid = len(recid)
        Nrecint = len(recidjmfint)
        
        recidrrm = np.intersect1d(recid, rrmprotsp)
        recjmfrrm = np.intersect1d(recjmf, rrmprotsp)
        recidjmfrrm = np.union1d(recidrrm, recjmfrrm)
        recidjmfintrrm = np.intersect1d(recidrrm, recjmfrrm)
        
        NrectotRRM = len(recidjmfrrm)
        NrecjmfRRM = len(recjmfrrm)
        NrecidRRM = len(recidrrm)
        NrecintRRM = len(recidjmfintrrm)
        
        
        recidkh = np.intersect1d(recid, khprotsp)
        recjmfkh = np.intersect1d(recjmf, khprotsp)
        recidjmfkh = np.union1d(recidkh, recjmfkh)
        recidjmfintkh = np.intersect1d(recidkh, recjmfkh)
        
        NrectotKH = len(recidjmfkh)
        NrecjmfKH = len(recjmfkh)
        NrecidKH = len(recidkh)
        NrecintKH = len(recidjmfintkh)
        
        rncmptjmf = np.unique(hfile[sfile[:, -1].astype(float) <= simcut,1])
        rncmptid = np.unique(ifile[ifile[:, -1].astype(float) >= idcut,1])
        
        motifjmf = np.unique(clusters[np.isin(clustnames,rncmptjmf)])
        motifid = np.unique(clusters[np.isin(clustnames,rncmptid)])
        
        Nrncmptjmf = len(rncmptjmf)
        Nrncmptid = len(rncmptid)
        Nmotifjmf = len(motifjmf)
        Nmotifid = len(motifid)
        
        wobj.write(specy.strip('/')+','+str(Nrbpprot)+','+str(NprotRRM)+','+str(NprotKH)+','+str(Nmeas)+','+str(NmeasRRM)+','+str(NmeasKH)+','+str(Nrectot)+','+str(Nrecjmf)+','+str(Nrecid)+','+str(Nrecint)+','+str(NrectotRRM)+','+str(NrecjmfRRM)+','+str(NrecidRRM)+','+str(NrecintRRM)+','+str(NrectotKH)+','+str(NrecjmfKH)+','+str(NrecidKH)+','+str(NrecintKH)+','+str(Nrncmptjmf)+','+str(Nrncmptid)+','+str(Nmotifjmf)+','+str(Nmotifid)+'\n')



