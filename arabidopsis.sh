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



# Use reconstructed motifs to infer RBPs that regulate RNA degradation rate
# Degradation rate was determined for tissues in Arabidopsis from intron and exon rpkm
# Degradation rate extracted by make_degradation-rate.py as the log(delta_exon)-log(delta_intron) - bias(delta_intro) 
# See for details: Alkallas, Rached, et al. "Inference of RNA decay rate from transcriptional profiling highlights the regulatory programs of Alzheimer’s disease." Nature communications 8.1 (2017): 1-11. https://www.nature.com/articles/s41467-017-00867-z

stability=${outdir}'StabilityAT'
mkdir $stability
degradationrate=${stability}'Degradation_rate.npz'
exonexpression=${stability}'Exon_expression_rate.npz'
atpwms=${stability}'Arabidopsis_thaliana_Reconstructionjplesvd122_svdsignificant_response122_confidence0.127.hmot'
at3putr=${stability}'Arabidopsis_thaliana.TAIR10.three_prime_utr.single_exon.longest_isoform.npz'
rbpexpression=${stability}'RBPs_Exon_expression_rate.npz'
atgeneset=${stability}'AT_highcorr_geneset.txt'
atnames=${stability}'Arabidopsis_thaliana_proteinnames.txt'
# Make RBP expression file
#python rbp_expression_file.py Arabidopsis_thaliana_Reconstructionjplesvd122_svdsignificant_response122_confidence0.127.hmot Arabidopsis_thaliana_proteinnames.txt Exon_expression_rate.npz


# Give tissue clusters names manually
tisclust='Tissues_custers.named.txt'
# Scan 3'UTR sequences with PWM for potential binding sites
python Scripts/scanFastaWithPWMregions.py $at2putr $atpwms AT.single_exon.longest_isoform_101pwmscan.npz --summax 0.1 --geneset $atgeneset
# Perform correlation analysis
python Scripts/GeneExpressioncorrelationAnalysis.py $degradationrate AT.single_exon.longest_isoform_101pwmscan-summaxgt0.1.npz $rbpexpression --correlationanalysis all 500 False optimal correlation --tissuemean SraRunTable.txt --lognorm_rbpexpression --outname stability --npz 
# Scan 100 shuffled 3'UTRs
for i in {0..100}; do submitjob python scanFastaWithPWMregions.py $at2putr $atpwms AT.single_exon.longest_isoform_101pwmscan.npz --summax 0.1 --geneset AT_highcorr_geneset.txt --random_shuffle $i; done
# Perform correlation analysis for shuffled sets
for i in {0..100}; do submitjob python Scripts/GeneExpressioncorrelationAnalysis.py $degradationrate AT.single_exon.longest_isoform_101pwmscan-summaxgt0.1_random_shuffle${i}.npz $rbpexpression --correlationanalysis all 500 False optimal correlation --tissuemean ../SraRunTable.txt --lognorm_rbpexpression --outname stability ; done
# Determine RBPs that signifianctly change the distribution
python Scripts/z_of_zscore.py AT.longest_isoform_101pwmscan_full-summaxgt0.1_random_shuffle:-out-tissuemean_correlation_all_optimal500.0False_correlationanalysis.txt 100 AT.single_exon.longest_isoform_101pwmscan_full-summaxgt0.1-out-tissuemean_correlation_all_optimal500.0False_correlationanalysis.txt 0.05 0.1
# Determine target sequences that correlate with RBPs expression (bona fide targets)
python Scripts/target_sets.py AT.single_exon.longest_isoform_101pwmscan_full-summaxgt0.1-out-tissuemean_correlation_all_optimal500.0False_correlationanalysis.npz AT.single_exon.longest_isoform_101pwmscan_full-summaxgt0.1-out-tissuemean_correlation_all_optimal500.0False_correlationanalysis-fdr0.1_zofz0.05.pvals 1,-1 $atnames $atpwms --cut 0.2

atsig='AT.single_exon.longest_isoform_101pwmscan_full-summaxgt0.1-out-tissuemean_correlation_all_optimal500.0False_correlationanalysis-fdr0.1_zofz0.05.pvals'
atout='AT.single_exon.longest_isoform_101pwmscan_full-summaxgt0.1-out-tissuemean_correlation_all_optimal500.0False_correlationanalysis.npz'



