# RBP sequence specificity
This folder contains the supplementary code for "Reconstructing the sequence preferences of RNA binding proteins across eukaryotes".

We used derived a linear embedding model to understand the relationship between protein sequence patterns and RNA sequence specificity. 

Dependencies:
python 2.7
scikit-learn, numpy, scipy, matplotlib, statsmodels, seaborn, biopython, etc. 
hmmer (http://hmmer.org/)

Before reconstructing the figures in the paper, please execute the data processing scripts in the following order:

1.A rncmpt_data.sh
1.B fig1.sh
2.A performance_calc.sh
2.B interface_importance.sh
2.C fig2.sh
3.A jple_reconstruction.sh
3.B fig3.sh
4.A cisbp-recstats.sh
4.B arabidopsis.sh
4.C fig4.sh
S.1 cisbp_reconstruction.sh

To reconstruct the individual figures in the paper run figX.sh.



