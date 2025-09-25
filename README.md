# Old code collection for: [A resource of RNA-binding protein motifs across eukaryotes reveals evolutionary dynamics and gene-regulatory function](https://www.nature.com/articles/s41587-025-02733-6)

## Please use [https://github.com/morrislab/jple](https://github.com/morrislab/jple/) instead

If you want to use JPLE to train your own model for predicting sequence specificity profiles from protein sequence, or infer RNA sequence specificity for RRM- or KH-domain using RNA-binding proteins, please refer to this repository: [https://github.com/morrislab/jple](https://github.com/morrislab/jple/)

All the measured and inferred RNA sequence specificities can be found at [cisbp.org](https://cisbp.org/)

### Content

Maintainance for this repository is limited. It contains collections of scripts and CLIs that were used for analyses during my PhD thesis ["Inferring RNA Sequence Specificities from Protein Sequences to Characterize Post-Transcriptional Regulation in Eukaryotes"](https://www.proquest.com/openview/d075b03685b79e4426cb11a18d9f19ad/1?pq-origsite=gscholar&cbl=18750&diss=y) 

Recommended: install anaconda and create virtual environment with python 2.7 adding dependencies listed in dependencies.txt

Additional requirements: 

- python3 (only agglomerative_clustering.py), scikit-learn==0.23.2)
- Hmmer (http://hmmer.org/)
- conservation_code (Capra JA and Singh M. Predicting functionally important residues from sequence conservation. Bioinformatics, 23(15):1875-82, 2007.
(https://compbio.cs.princeton.edu/conservation/))
- pymol (https://pymol.org/2/) to visualize individual pdbs

### Reproduce thesis results:

Before reconstructing the figures with (fig1-5.sh), execute data processing scripts (rncmpt_data.sh, performance_calc.sh, interface_importance.sh, jple_reconstruction.sh, cisbp-recstats.sh, arabidopsis.sh)

The following order is required:

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

### License

All the code is licensed under the BSD 3-Clause License. See [LICENSE](https://github.com/lxsasse/RBPbinding/blob/main/LICENSE.md) for details.



