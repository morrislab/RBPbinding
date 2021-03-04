This directory contains python scripts that were used to generate the data and figures

# scriptname: brief description

agglomerative_clustering.py : clusters proteins on binding similarity, sequence identity or sub-clusters cluster from other metric (python 3, because scikit-learn changed)

Eigenvalue.py : Visualizes similarities and differences between SVD on joint and individual protein k-mers and binding enrichments

GeneExpressioncorrelationAnalysis.py : Detect differential correlation between bound and unbound sets of genes

kmer_map2seq.py : Align protein k-mers back to protein sequence and sum up k-mer scores for interface residue scores

protein-features_beta.py : Generate protein k-mer profile for given protein sequences

scanFastaWithPWMregions.py : Scan 3'UTR sequences with binding PSSMs

Specificity_predictors.py : Predict RNA sequence specficity with various methods (JPLE, LR, KNN, Affinity regression)

top_identity_to_trainset.py : Predict RNA sequence specificity from most similar protein sequence
