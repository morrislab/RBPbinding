import numpy as np
import sys, os

dataarray = np.genfromtxt(sys.argv[1], delimiter = ';', dtype = str)
species = open(sys.argv[1], 'r').readline().strip().split(';')[16:]
# Make file for species specific values:
# Number of each remaining RSSGs orgin
# Number of each remaining evolutionary types
# Number for KH and RRM seperately
cluster = dataarray[:, 16:].astype(int)

Existing_clades, exindex = np.unique(dataarray[:,14], return_index = True)
Ageclade = dataarray[exindex,8]
Existing_clades = list(Existing_clades[np.argsort(-Ageclade.astype(float))])

Existing_evolution = list(np.unique(dataarray[:,12]))
EVnames = ['Leca_Ancient', 'Duplication', 'Parental_Specialization', 'Sibling_Specialization', 'Insecure']

Domain_class = dataarray[:,4]
khspecialfrac = float(np.sum(dataarray[Domain_class == 'KH',12] == '2'))/float(np.sum(np.isin(dataarray[Domain_class == 'KH', 12], ['1','2'])))
rmmspecialfrac = float(np.sum(dataarray[Domain_class == 'RRM',12] == '2'))/float(np.sum(np.isin(dataarray[Domain_class == 'RRM', 12], ['1','2'])))
print 'Fraction of speciation for RRM and KH domains', khspecialfrac,rmmspecialfrac 
print 'Number of RRM and KH RSSGs', int(np.sum(Domain_class == 'RRM')), int(np.sum(Domain_class == 'KH'))


#for d, dom in enumerate(Domain_class):
    #Domain_class[d] = dom.split('-')[0]

specoutmatrix = np.zeros((len(species), len(Existing_clades)+len(EVnames)*3), dtype = int)
for s, spec in enumerate(species):
    chosencl = cluster[:,s]>=1
    clclade, clcladeN = np.unique(dataarray[chosencl, 14], return_counts = True)
    clevo, clevoN = np.unique(dataarray[chosencl, 12], return_counts = True)
    for c, cl in enumerate(clclade):
        specoutmatrix[s, Existing_clades.index(cl)] = clcladeN[c]
    for c, cl in enumerate(clevo):
        specoutmatrix[s, len(Existing_clades) + Existing_evolution.index(cl)] = clevoN[c]
    for d, dc in enumerate(['RRM', 'KH']):
        chosendcl = chosencl * (Domain_class == dc)
        clevo, clevoN = np.unique(dataarray[chosendcl, 12], return_counts = True)
        for c, cl in enumerate(clevo):
            specoutmatrix[s, len(Existing_clades) + len(Existing_evolution)*(1+d)+Existing_evolution.index(cl)] = clevoN[c]
Header = np.concatenate([['Species'], Existing_clades, EVnames, ['RRM_'+EVnames[i] for i in range(5)], ['KH_'+EVnames[i] for i in range(5)]])
np.savetxt(os.path.splitext(sys.argv[1])[0]+'_species_numbers.txt', np.append(np.array(species)[:,None], specoutmatrix.astype(int),axis = 1), header = ';'.join(Header), delimiter = ';', fmt = '%s')

            
            
            
            
            
            
            
                  
            
            
