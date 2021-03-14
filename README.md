# RBP sequence specificity
This directory contains the code for "Reconstructing sequence specificities of RNA binding proteins across eukaryotes".

We developed a linear embedding model to understand the relationship between protein sequence patterns and RNA sequence specificity. 

Recommended: create virtual python 2.7 environment and install dependencies in dependencies.txt

Additional requirements: 

- python3 (only agglomerative_clustering.py), scikit-learn==0.23.2)
- Hmmer (http://hmmer.org/)
- conservation_code (Capra JA and Singh M. Predicting functionally important residues from sequence conservation. Bioinformatics, 23(15):1875-82, 2007.
(https://compbio.cs.princeton.edu/conservation/))
- pymol (https://pymol.org/2/) to visualize individual pdbs

Note: 
To execute "full" pipeline (i.e. every intermediate step), very long running times and storage is required. Parallel execution recommended!
To run modify intermediate results, set $full=1 in bash scripts.

RUN:
Before reconstructing the figures with (fig1-5.sh), execute data processing scripts (rncmpt_data.sh, performance_calc.sh, interface_importance.sh, jple_reconstruction.sh, cisbp-recstats.sh, arabidopsis.sh)

The following order is recommended/required:

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

cisbp_reconstruction.sh executes scripts to locally generate data available on http://cisbp-rna.ccbr.utoronto.ca/, e.g. PWMs jpgs for confidentily reconstructed specificities. 




