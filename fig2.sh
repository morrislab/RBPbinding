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


Imask=${outdir}'RRM_Imask/'

jplepredpdb=${Imask}'Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequenceRncmpt_templates-id50.0-msim3.5_all.txt'
conspredpdb=${Imask}'Conservation_RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated.domain.homolog30_homologs_maxsim0.95_minhom0.3Rncmpt_templates-id50.0-msim3.5_all.txt'
rfpredpdb=${Imask}'rf_complete_wz5sum-Nonexnorm-testinterfaceRncmpt_templates-id50.0-msim3.5_all.txt'


##### Figure 2 plots ######
# Figure 2B, grey UMAP embedding
python Scripts/motifset_cluster_4.py --plotfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz 1 umap 2 --reduceparams 14,cosine,0.7,1.,42 --plotevo --figuresettings 6 6 0.3 100 False 0 0 0 --savefig ${fig2}JPLElatent_umap.jpg
# Figure 2C, Pearson recall curve for sequence identity and JPLE
python Scripts/precision_recall_predictions.py 2 --scores -${outdir}Performance/JPLE/jplesvd122_set_complete__svdsignificant_response122.0_maplsq_declocal_-testset_latentstats.dat,-1 ${outdir}Performance/SeqID/RNCMPT_unique_experiments_RRM,KH_testsetSingle_Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id_bid.dat,-1 --realval ${outdir}Performance/JPLE/jplesvd122_set_complete__svdsignificant_response122.0_maplsq_declocal_-testset_profile_pcorrelation.dat,-2 ${outdir}Performance/SeqID/RNCMPT_unique_experiments_RRM,KH_testsetSingle_Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id_bid-Zscores_420_origrncmpt_pearson.dat,-1 --setnames "JPLE confidence, SeqID" --legend --savefig ${fig2}Pearson_recall_jple-seqId.jpg --definecolors 1,0 --ycut 0.7
# Figure 2D, AUROC curve plot
python Scripts/kmer_AUC.2.py ${Imask}Rncmpt_templates-id50.0-msim3.5_all.txt $jplepredpdb $conspredpdb JPLE,Conservation --boxroc --boxpr --diffscatter --combine --savefig ${fig2}AUC_pdb_jple_vs_conservation.jpg --saveresults
# Figure 2D, structure plot
echo 'Save pymol figure as png'
pymol ${Imask}RNCMPT00121_to_4ed5_B_CRncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.pdb ${Imask}showpdb.pml B C


##### Supplementary Figures Figure2
# Plot k-mer overlap between KH and RRM
python Scripts/venn_kmeroverlap.py ${indir}${masterfile} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz RRM,KH --outdir $supfig
# Estimate noise in specificty measurements and determine optimal number of eigenvectors
python Scripts/Eigenvalue.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz ${indir}${rrmkhlist} 0.95 ${supfig}
# Comparison to other methods
# Pearson-recall to Protein vector and Affinity regression
python Scripts/precision_recall_predictions.py 4 --scores ${outdir}Performance/SeqID/RNCMPT_unique_experiments_RRM,KH_testsetSingle_Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id_bid.dat,-1 -${outdir}Performance/JPLE/jplesvd122_set_complete__svdsignificant_response122.0_maplsq_declocal_-testset_latentstats.dat,-1 ${outdir}Performance/AR/ar_set_complete_AR-P0.9_YD0.95_centFalse_alph0.1_ridgeFalse_D5-directfit_highcorrelation.dat,-1 -${outdir}Performance/NN/nn1cosine_set_complete__KNN1.0-cosine-equal-testset_latentstats.dat,-1 --realval ${outdir}Performance/SeqID/RNCMPT_unique_experiments_RRM,KH_testsetSingle_Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id_bid-Zscores_420_origrncmpt_pearson.dat,-1  ${outdir}Performance/JPLE/jplesvd122_set_complete__svdsignificant_response122.0_maplsq_declocal_-testset_profile_pcorrelation.dat,-2 ${outdir}Performance/AR/ar_set_complete_AR-P0.9_YD0.95_centFalse_alph0.1_ridgeFalse_D5-testset_profile_pcorrelation.dat,-2 ${outdir}Performance/NN/nn1cosine_set_complete__KNN1.0-cosine-equal-testset_profile_pcorrelation.dat,-2 --setnames "SeqID, JPLE confidence, AR closest specificity, Nearest neighbor 5mer" --legend --savefig ${supfig}Pearson_recall_jple-AR-nn5mer-seqId.jpg --definecolors 5,1,3,2
# Scatter-plot to DeepL features
python Scripts/scatter_predictions.py --files ${outdir}Performance/DeepL/lrunirep_set_complete__LR0.0001fiTrue-testset_profile_pcorrelation.dat ${outdir}Performance/JPLE/jplesvd122_set_complete__svdsignificant_response122.0_maplsq_declocal_-testset_profile_pcorrelation.dat --definecolums -2 -2 --minlength 354 --setnames UnirepDLfeatures JPLE  --plotlim 0,1.025 --label 'Pearson' --savefig ${supfig}Scatter_Pearson_JPLE-unirep.jpg
python Scripts/scatter_predictions.py --files ${outdir}Performance/DeepL/lrbertrep_set_complete__LR0.0001fiTrue-testset_profile_pcorrelation.dat ${outdir}Performance/JPLE/jplesvd122_set_complete__svdsignificant_response122.0_maplsq_declocal_-testset_profile_pcorrelation.dat --definecolums -2 -2 --minlength 354 --setnames UnirepDLfeatures JPLE  --plotlim 0,1.025 --label 'Pearson' --savefig ${supfig}Scatter_Pearson_JPLE-bert.jpg
# Individual Area under the receiver operating curves against conservation, and interface visualizations
python Scripts/kmer_AUC.2.py  ${Imask}Rncmpt_templates-id50.0-msim3.5_all.txt $jplepredpdb $conspredpdb JPLE,Conservation  --prplot --rocplots --interface_visualization 5. --combine --savefig ${supfig}AUC_pdb_jple_vs_conservation.jpg
# Scatter plot comparison to SVM trained on structures
python Scripts/kmer_AUC.2.py ${Imask}Rncmpt_templates-id50.0-msim3.5_all.txt $jplepredpdb $rfpredpdb JPLE,Conservation --boxroc --boxpr --diffscatter --combine --savefig ${supfig}AUC_pdb_jple_vs_rf.jpg



