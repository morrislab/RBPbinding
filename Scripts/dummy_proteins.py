import numpy as np
import sys, os

species = np.genfromtxt(sys.argv[-1], dtype = str)
sidcut = float(sys.argv[-2])

kmers = []
features = []
names = []
for s, spec in enumerate(species):
    sidfile = np.load(spec+'/'+spec+sys.argv[1])
    maxid = np.amax(sidfile['identmat'], axis =1)
    keepprot = sidfile['names2']
    keepprot = keepprot[maxid <= sidcut]
    for k, kep in enumerate(keepprot):
        keepprot[k] = kep.split('__')[0]
    sfefile = np.load(spec+'/'+spec+sys.argv[2])
    protnames = sfefile['protnames']
    pmask = np.isin(protnames, keepprot)
    if np.sum(pmask) > 0:
        kmers.append(sfefile['kmers'])
        names.append(protnames[pmask])
        pfeat = sfefile['features']
        features.append(pfeat[:, pmask])

nkmers = np.unique(np.concatenate(kmers))
names = np.concatenate(names)
nfeatures = []
for f, feat in enumerate(features):
    nfeat = np.zeros((len(nkmers),len(feat[0])))
    nfeat[np.isin(nkmers, kmers[f])] = feat
    nfeatures.append(nfeat)

nfeatures = np.concatenate(nfeatures, axis = 1)
print np.shape(nfeatures), len(names), len(nkmers)
np.savez_compressed('Dummy_RRMKH_49species_5pm1_features_idmax'+str(sidcut)+'.npz', features = nfeatures, expnames = names, protnames = names, kmers = nkmers)


