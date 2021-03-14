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

seqid=${outdir}'Rncmpt.aaseq.ext15_combined_lnormmin_alignlocal_id.txt'
proteinnames=${indir}'RNCMPT_protein_names.txt'

#### Figure 3 ###
# 3A Umap of latent space with sequence identity connections
python Scripts/motifset_cluster_4.py --plotfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz 1 umap 2 --reduceparams 14,cosine,0.7,1.,42 --clusterassignment centrecolor ${outdir}Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6.npz --largecluster 5 --color_motifcomposition ${outdir}Zscores_420_origrncmpt_pwm.hmot pwm manual --annotate coloring outside 5 --motifshape ${indir}Rncmpt.aaseq_domainclass.txt --connect $seqid 40 --plotevo --savefig ${fig3}JPLElatent_with_sequenceidentity_motifcluster5.jpg
# Returns UMAP and PWMs for clusters, plot clust pwms with 
python Scripts/pylogo.py ${outdir}Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6-motifsetvisual-minsize5-pwms.txt --infocont --removeaxis --transparent --format png

# 3B Umap of latent space for eigenvector 1 and 2 adn 3
python Scripts/motifset_cluster_4.py --plotfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz 1 umap 2 --reduceparams 14,cosine,0.7,1.,42 --color_factorplot 2 None True True $proteinnames --savefig ${fig3}JPLElatent_umapeigenvectors.jpg

# 3C Residue importance profiles for SART3, sart3(Tn), Rnp4F(Dm)
python Scripts/interface_secstruc_map.py ${outdir}RRM_Imask/Rncmpt.RRMKH.ss.split.rn.fasta ${outdir}Rncmpt.aaseq.ext15_domain_fusedjplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta --fullseqloc ${indir}Rncmpt.aaseq_pseq.fasta --ylimmean --savefig ${fig3}Rncmpt.aaseq.ext15_domain_fusedjplesvd122_PrecZ.jpg --proteinset RNCMPT00060,RNCMPT00064,RNCMPT00774
python Scripts/interface_importance_comparison.py ${outdir}Rncmpt.aaseq.ext15_domain_fusedjplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta RNCMPT00774,RNCMPT00064,RNCMPT00060 --yzscale --savefig Figures/Figure3/Sart3comparison.jpg

# 3D Distribution similar and dissimilar latent distances for JPLE and sequence identity 
python Scripts/precision_recall_similarities.py 1 --scorenpz 354 training -${outdir}Performance/JPLE/jplesvd122_set:_svdsignificant_response122.0_maplsq_declocal_-testset_latentdists.npz:cosine --Identities $seqid --realval ${outdir}Zscores_420_origrncmpt_pearson.mat --setnames "JPLE latent distance,Sequence identity" --cutoff 0.6 --savefig ${fig3}Similaritiescut0.6_lve1out-trainset_cosine_JPLE-ID.jpg --definecolors 0,1 --densities 0.6 0.6

# 3E UMAP embedding with proteins from 49 species
python Scripts/motifset_cluster_4.py --plotfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz 1 umap 2 --reduceparams 14,cosine,0.7,1.,42 --clusterassignment centrecolor ${outdir}Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6.npz --largecluster 5 --color_motifcomposition ${outdir}Zscores_420_origrncmpt_pwm.hmot pwm manual --annotate coloring outside 5 --plotevo --secondfeatures ${outdir}JPLE_RRMKHcomplete/RecomputeDummyjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdecglobal-testset_Platent.npz  --savefig ${fig3}JPLElatent_withDummy_motifcluster5.jpg
# 3F Finding correct sequence for unknown binding specificity
python Scripts/sequence_rank.py ${outdir}Performance/JPLEunknown/jplesvd122_set_complete_svdsignificant_response122.0_maplsq_declocal_Dummy_RRMKH_49species_5pm1_features_rank.txt,${outdir}Performance/LRunknown/lr_set_complete_LR0.0fiFalseDummy_RRMKH_49species_5pm1_features_rank.txt 4050 'JPLE latent,LR reverse' --savefig ${fig3}JPLE_LR_sequence_search.jpg

#### Supplementary Figures Figure 3 
# Umap with eigenvector coloring of 30 eigenvectors
python Scripts/motifset_cluster_4.py --plotfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz 1 umap 2 --reduceparams 14,cosine,0.7,1.,42 --color_factorplot 30 None False True $proteinnames --savefig ${supfig}JPLElatent_umapeigenvectors.jpg
# Umap with location of Rnp4F, sart3, SART3
python Scripts/motifset_cluster_4.py --plotfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz 1 umap 2 --reduceparams 14,cosine,0.7,1.,42 --clusterassignment centrecolor ${outdir}Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6.npz --largecluster 5 --color_motifcomposition ${outdir}Zscores_420_origrncmpt_pwm.hmot pwm manual --annotate coloring outside 5 --motifshape ${indir}Rncmpt.aaseq_domainclass.txt --connect $seqid 40 --plotevo --savefig ${fig3}JPLElatent_with_sequenceidentity_motifcluster5_sart3.jpg --locate_protein $proteinnames RNCMPT00060,RNCMPT00064,RNCMPT00774
# Eigenvector distribution of entries and summary of content
python Scripts/Factor-visualization-beta.py ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef.npz 0.96 0.44 30 --makelist --bothends --zkmerchoice 10 6 5 --pkmerchoice 20 4 3 --infocont --removespines --plot_scoredistribution --savefig --outdir ${supfig}

