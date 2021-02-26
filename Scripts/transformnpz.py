improt numpy as np

f = np.load('Zscores_420_origrncmpt.npz')
zscores = f['zscores']
kmers = f['kmers']
expnames = f['expnames']

np.savetxt('Zscores_420_origrncmpt.txt', np.append(kmers.reshape(-1,1), zscores, axis = 1).astype(str), header = ' '.join(expnames), fmt = '%s')




