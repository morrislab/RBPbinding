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



#### Pre-calculations for Figure 2, 3, 4, 5####
# compute k-mer count
# gapped 5-mers for JPLE
python Scripts/protein-features_beta.py --mprocessing 4 --domainfasta ${outdir}Rncmpt.aaseq.ext15_domain_fused.fasta --k 5 --outname ${outdir}Rncmpt.aaseq.ext15_domain_fused_5pm1 --posmismatchkmer 1
# ungapped 4mers for Affinity regression
python Scripts/protein-features_beta.py --mprocessing 4 --domainfasta ${outdir}Rncmpt.aaseq.ext15_domain_fused.fasta --k 4 --outname ${outdir}Rncmpt.aaseq.ext15_domain_fused
## gapped 5-mers with homologs 
cisbpext15='CisBP_RRMKH.ext15_domain_concatenated.fasta'
# use phmmer to scan CisBP for similar sequences
phmmer -o Rncmpt.aaseq.ext15_domain_fused.homologs.out -E 1e-15 --domE 1e-15 --incE 1e-15 Rncmpt.aaseq.ext15_domain_fused.fasta $cisbpext15
# extract sequences with minimum and maximum sequence identity to RNCMPT measurements
python Scripts/Hmmerout_getorthologs.py Rncmpt.aaseq.ext15_domain_fused.homologs.out Rncmpt.aaseq.ext15_domain_fused.fasta 0.99 0.5 Rncmpt.aaseq.ext15_domain_fused
# extract gapped 5-mers from original sequences and homologs sequences to train semi-supervised model
python Scripts/protein-features_beta.py --mprocessing 8 --domainfasta Rncmpt.aaseq.ext15_domain_fused.fasta --domain_orthologs Rncmpt.aaseq.ext15_RRMKH_fused_homologs_maxsim0.99_minhom0.5.fasta --k 5 --outname Rncmpt.aaseq.ext15_domain_fused_5hom99-50 --posmismatch 1


# Define leave-one-out training and test sets for RRM and KH domains to assess selectivity and sensitivity 
python Scripts/reconstruction_trainingsets.py ${indir}${rrmkhlist} None None None None False Single
# define leave-one-out training and test sets for distinct/unknown specificities to assess desired sequence search
python Scripts/reconstruction_trainingsets.py ${indir}${rrmkhlist} ${outdir}Zscores_420_origrncmpt_pearson.mat 0.6 None None False Single

# Run JPLE on full data set for UMAP visualization and RBR interface predictions
compdir='JPLE_RRMKHcomplete'
mkdir $compdir
# Supervised version
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5pm1_5mer_features.npz --proteinset ${outdir}${rrmkhlist} --JPLE svd 122 significant_response lsq global --savetestcorrelation --savetopintersection --savetestprofileslight --savelatentstats --savemodel --savetestprofiles --savelatenttest --save_reconstructionP --normP2 --normY2 --outname ${compdir}jplesvd122
# Semi-supervised version
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5hom99-50_5mer_features.npz --proteinset ${outdir}${rrmkhlist} --JPLE svd 0.95 significant_response lsq global --savetestcorrelation --savetopintersection --savetestprofileslight --savelatentstats --savemodel --savetestprofiles --savelatenttest --save_reconstructionP --normP2 --normY2 --outname ${compdir}jplesvd0.95hom


## Comparison of performance of predictors
perform='Performance/'
mkdir $perform
mkdir ${perform}JPLE/
mkdir ${perform}AR/
mkdir ${perform}NN/
mkdir ${perform}DeepL/
mkdir ${perform}JPLE5hom50/

for i in {0..354}
do
# Run JPLE for all leave-one out sets
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5pm1_5mer_features.npz --trainingset $i ${outdir}${rrmkhlist%.list}_testsetSingle.list --JPLE svd 122 significant_response lsq local --savetestcorrelation --savetopintersection --saveclosestcorrelation --savetopcorrelation --savetestprofileslight --savelatentstats --normP2 --normY2 --outname ${perform}JPLE/jplesvd122_set${i}

# Run Semi-supervised JPLE for leave-one out
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir} Rncmpt.aaseq.ext15_domain_fused_5pm1hom50-99_5mer_features.npz --trainingset $i ${outdir}${rrmkhlist%.list}_testsetSingle.list --JPLE svd 122 significant_response lsq local --savetestcorrelation --savetopintersection --saveclosestcorrelation --savetopcorrelation --savetestprofileslight --savelatentstats --savemodel --normP2 --normY2 --outname Performance/JPLE5hom50/jpletrainset${i}
# Assess semi-supervised trained model with non-homology features
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5pm1hom50-99_5mer_features.npz --JPLE svd 122 significant_response lsq local --savetestcorrelation --savetopintersection --saveclosestcorrelation --savetopcorrelation --savetestprofileslight --savelatentstats --normP2 --normY2 --recompute Performance/JPLE5hom50/jpletrainset${i}_svdsignificant_response122.0_maplsq_declocal__coef.npz Rncmpt.aaseq.ext15_domain_fused_5pm1_5mer_features.npz Zscores_420_origrncmpt.txt testprots Performance/JPLE5hom50/Recompute

# Run Affinity regression for all leave-one out sets
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_4mer_features.npz --trainingset $i ${outdir}${rrmkhlist%.list}_testsetSingle.list --Affinity_regression 5 0.9 0.95 False 0.1 False --savetestcorrelation --savetopintersection --saveclosestcorrelation --savetopcorrelation --save_predicted_similarity --savetestprofileslight --normP2 --normY2 --outname ${perform}AR/ar_set${i}

# Run Nearest Neighbor for all leave-one out sets
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz --trainingset $i ${outdir}${rrmkhlist%.list}_testsetSingle.list --knearestneighbor 1 cosine equal --savetestcorrelation --savetopintersection --saveclosestcorrelation --savetopcorrelation --savetestprofileslight --savelatentstats --outname ${perform}NN/knn1cosine_set${i}

# Run Regression with Deep learning features for all leave one out sets
python Scripts/Specificity_predictors.py ${indir}${zscores} ${indir}TAPE_unirep.fasta.avgsum.npz --trainingset $i ${outdir}${rrmkhlist%.list}_testsetSingle.list --LinearRegression 0.0001 True False --savetestcorrelation --savetopintersection --saveclosestcorrelation --savetopcorrelation --savetestprofileslight --outname ${perform}DeepL/lrunirep_set${i}
python Scripts/Specificity_predictors.py ${indir}${zscores} ${indir}TAPE_bert.fasta.avgsum.npz --trainingset $i ${outdir}${rrmkhlist%.list}_testsetSingle.list --LinearRegression 0.0001 True False --savetestcorrelation --savetopintersection --saveclosestcorrelation --savetopcorrelation --savetestprofileslight --outname ${perform}DeepL/lrbert_set{i}
done

# concatenate single output files for analysis
python Scripts/concatenatefiles.py ${perform}JPLE/jplesvd122 _set .dat complete --minnum 355
python Scripts/concatenatefiles.py ${perform}JPLE5hom50/Recompute _set .dat complete --minnum 355
python Scripts/concatenatefiles.py ${perform}AR/ar _set .dat complete --minnum 355
python Scripts/concatenatefiles.py ${perform}NN/knn1cosine _set .dat complete --minnum 355
python Scripts/concatenatefiles.py ${perform}DeepL/lrunirep _set .dat complete --minnum 355
python Scripts/concatenatefiles.py ${perform}DeepL/lrbert _set .dat complete --minnum 355

# Run sequence identity for all leave-one out sets as comparison
mkdir ${perform}SeqID
python Scripts/top_identity_to_trainset.py 0-355 RNCMPT_unique_experiments_RRM,KH_testsetSingle.list Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt --similaritymat ${outdir}${zscores%.txt}_pearson.mat --savetxt ${perform}SeqID/


