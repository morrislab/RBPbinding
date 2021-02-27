# RBP sequence specificity
This folder contains the supplementary code for "Reconstructing the sequence preferences of RNA binding proteins across eukaryotes".

We used derived a linear embedding model to understand the relationship between protein sequence patterns and RNA sequence specificity. 

Dependencies:

python 2.7

scikit-learn, numpy, scipy, matplotlib, statsmodels, seaborn, biopython, etc.

Hmmer (http://hmmer.org/)

Capra JA and Singh M. Predicting functionally important residues from sequence conservation. Bioinformatics, 23(15):1875-82, 2007.
(https://compbio.cs.princeton.edu/conservation/)

Before reconstructing the figures in the paper, please execute the data processing scripts in the following order:

1. rncmpt_data.sh
2. fig1.sh
3. performance_calc.sh
4. interface_importance.sh
5. fig2.sh
6. jple_reconstruction.sh
7. fig3.sh
8. cisbp-recstats.sh
9. arabidopsis.sh
10. fig4.sh
11. cisbp_reconstruction.sh

To reconstruct the individual figures in the paper run figX.sh.



