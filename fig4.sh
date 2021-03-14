masterfile='S1.In_Progress.v6.For_all.DR-2.csv'
# transform Zscore npz into txt file to read and read in 
zscores='Zscores_420_origrncmpt.txt'


rrmkhlist='RNCMPT_unique_experiments_RRM,KH.list'
khlist='RNCMPT_unique_experiments_KH.list'
rrmlist='RNCMPT_unique_experiments_RRM.list'

highestzscoreset='RNCMPT_unique_experiments_zcoresmax381.list'

proteinnames='RNCMPT_protein_names.txt'
domainclass='Rncmpt.aaseq_domainclass.txt'


indir='Data/'
outdir='Outdir/'

fig='Figures/'

fig1=${fig}'Figure1/'
fig2=${fig}'Figure2/'
fig3=${fig}'Figure3/'
fig4=${fig}'Figure4/'
fig5=${fig}'Figure5/'

supfig=${fig}'Supplementary/'

stability=${outdir}'Arabidopsis_thaliana/'
recdir=${outdir}'Reconstruction/'

cisbprecstats=${recdir}'CisBP_reconstruction-stats-jple0.127-id70.0-cisbp2-Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6_nspecies691.csv'
cisbprecsoft=${recdir}'CisBP_reconstruction-stats-jple0.127-id40.0-cisbp2-Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6_nspecies691.csv'
cisbprecstats06=${recdir}'CisBP_reconstruction-stats-jple0.127-id70.0-cisbp06-Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6_nspecies691.csv'

atsig=${stability}'AT.single_exon.longest_isoform_101pwmscan-summaxgt0.1-stability-tissuemean_correlation_all_optimal500.0False_correlationanalysis-fdr0.1_zofz0.05.pvals'
atout=${stability}'AT.single_exon.longest_isoform_101pwmscan-summaxgt0.1-stability-tissuemean_correlation_all_optimal500.0False_correlationanalysis.npz'
atpwm=${stability}'Arabidopsis_thaliana_Reconstructionjplesvd122_svdsignificant_response122_confidence0.127.hmot'
attargets=${stability}'AT.single_exon.longest_isoform_101pwmscan-summaxgt0.1-stability-tissuemean_correlation_all_optimal500.0False_correlationanalysis-fdr0.1_zofz0.05_cut0.2-targetgenes.npz '
atnames=${stability}'Arabidopsis_thaliana_proteinnames.txt'

#### Figure 4 #####
# Fig 4A
python Scripts/compare_version-coverage-2.py $cisbprecstats $cisbprecstats06 --specieslist ${recdir}species_handcraft_measured2.nwk.list --comparetoid --legend --savefig ${fig4}Coverage_Cisbpversion-leg.jpg
# Fig 4B, Coverage JPLE versus ID
python Scripts/coverage-3.py $cisbprecstats --specieslist ${recdir}specieslist.dat --kingdoms ${recdir}speciesclades.txt fraction --figsize 2.3,3 --savefig ${fig4}Coverage_JPLEvsID70.jpg --legend
# Fig 4C, Coverage JPLE versus ID scatterplot
python Scripts/coverage-3.py  $cisbprecstats --specieslist ${recdir}species_handcraft_measured2.nwk.list --kingcolor ${recdir}species_handcraft_measured2.nwk.txt 1 --figsize 3,3 --savefig ${fig4}Coverage_JPLEvsID70_49species_scatter.jpg --scattercomparison
# Fig 4D
python Scripts/plotzscore.py $atsig 1,2,-1 ${stability}Arabidopsis_thaliana_proteinnames.txt --savefig --outdir ${fig4}
python Scripts/plotzscore.py $atsig 1,2,-1 ${stability}Arabidopsis_thaliana_proteinnames.txt --savefig --legend --outdir $supfig
python Scripts/pylogo.py $atpwm --infocont --removeaxis --transparent --format png
# Fig 4E
python Scripts/rbp-expression.py ${stability}RBPs_Exon_expression_rate.npz ${stability}Degradation_rate.npz $atsig 1,-1 $atnames ${atpwm} ${stability}SraRunTable.txt --clustertissue ${stability}Tissues_custers.named.txt  --fullsimilarity --outdir ${fig4} --savefig
# Fig 4F
python Scripts/goterm-targets.py $attargets $atpwm ${stability}ATH_GO_GOSLIM.npz --savefig --outdir ${fig4}

### Supplementary figures
# Reconstruction comparison between JPLE and ID for 49 species detailed
python Scripts/coverage-3.py $cisbprecstats --specieslist ${recdir}species_handcraft_measured2.nwk.list --figsize 3,13 --savefig ${supfig}Coverage_JPLEvsID70_49species.jpg --legend
# Scatter comparison for RRM and KH only
python Scripts/coverage-3.py $cisbprecstats --specieslist ${recdir}species_handcraft_measured2.nwk.list --kingcolor ${recdir}species_handcraft_measured2.nwk.txt 1 --figsize 3,3 --savefig ${supfig}Coverage_JPLEvsID70_49species_rrm.jpg --scattercomparison --rrm
# Scatter comparison for RRM and KH to 40% seq id cutoff
python Scripts/coverage-3.py $cisbprecsoft --specieslist ${recdir}species_handcraft_measured2.nwk.list --kingcolor ${recdir}species_handcraft_measured2.nwk.txt 1 --figsize 3,3 --savefig ${supfig}Coverage_JPLEvsID40_49species_rrm.jpg --scattercomparison --rrm

# Relationship between sequence identity and latent distance of all eukayotic RBPs
python Scripts/resconstruction-relations.2.py ${recdir}CisBP_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat ${recdir}CisBP_RRMKH_domain_combined_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq_maxid.txt 0.127 70 --savefig --outdir $supfig

# Relationship between sequence identity and latent distance for each training object
python Scripts/resconstruction-relations_individual.py ${recdir}CisBP_highid_jplemin_species.stats ${indir}RNCMPT_protein_names.txt ${outdir}Zscores_420_origrncmpt_pwm.motif ${indir}RBP_420domainclass.txt 0.127 70 --savefig --outdir $supfig
# Check whether degradation rates are reproducible across genes and maintain tissue similarity

python Scripts/GeneExpressioncorrelation.py ${stability}Degradation_rate.npz ${stability}SraRunTable.txt 'Stability' --outdir $supfig > ${stability}Tissue_clusters.txt
# Correlation distribution for bound and unbound sequences
python Scripts/binding-effect.py $atout $atsig 1,-1 $atnames  $atpwm --savefig --outdir $supfig



