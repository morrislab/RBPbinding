import numpy as np
import sys, os

f = np.load(sys.argv[1])
zscores = f['zscores']
kmers = f['kmers']
expnames = f['expnames']

np.savetxt(os.path.splitext(sys.argv[1])[0]+'.txt', np.append(kmers.reshape(-1,1), zscores, axis = 1).astype(str), header = ' '.join(expnames), fmt = '%s')




