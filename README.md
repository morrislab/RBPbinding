# RBP sequence specificity
This directory contains the code for "Reconstructing sequence specificities of RNA binding proteins across eukaryotes".

We developed a linear embedding model to understand the relationship between protein sequence patterns and RNA sequence specificity. 

Dependencies:

python 2.7 (mostly), python 3 (only once)

scikit-learn, numpy, scipy, matplotlib, statsmodels 0.10.0, seaborn, biopython, logomaker etc.

To execute "full" pipeline (i.e. every intermediate step, very long running time required, parallel execution recommended), 
also download and install:

Hmmer (http://hmmer.org/)

Capra JA and Singh M. Predicting functionally important residues from sequence conservation. Bioinformatics, 23(15):1875-82, 2007.
(https://compbio.cs.princeton.edu/conservation/)

pymol (https://pymol.org/2/)

Before reconstructing the figures with (fig1-5.sh), execute the data processing scripts (rncmpt_data.sh, performance_calc.sh, interface_importance.sh, jple_reconstruction.sh, cisbp-recstats.sh, arabidopsis.sh)

The following order is recommended

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
11. fig5.sh

cisbp_reconstruction.sh executes scripts to generate the data available on http://cisbp-rna.ccbr.utoronto.ca/




