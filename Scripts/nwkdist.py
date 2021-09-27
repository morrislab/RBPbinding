import numpy as np
import sys, os
from ete3 import ClusterTree

newick = open(sys.argv[1], 'r').readline().strip()

#tree = ClusterTree(newick)
#cophenetic_matrix, newick_labels = tree.cophenetic_matrix()
#cophenetic_matrix = pd.DataFrame(cophenetic_matrix, columns=newick_labels, index=newick_labels)
import dendropy
from dendropy import treecalc

tree = dendropy.Tree.get_from_path(sys.argv[1], "newick")
pdm = tree.phylogenetic_distance_matrix()

names = tree.taxon_namespace
distance = np.zeros((len(names), len(names)))


for i, t1 in enumerate(names):
    for j in range(i+1, len(names)):
        t2 = names[j]
        distance[i,j] = distance[j,i] = pdm(t1,t2)/2.
        print("Distance between '%s' and '%s': %s" % (t1, t2, distance[i,j]))

names = np.array(names, dtype = str)
for n, name in enumerate(names):
    names[n] = name.replace(' ', '_').strip("'")

print( '\n'.join(np.array(names, dtype = str)))

print('\nSAVED as', sys.argv[1]+'distance_matrix.txt\n')
np.savetxt(sys.argv[1]+'distance_matrix.txt', distance, header = ' '.join(np.array(names, dtype = str)))



