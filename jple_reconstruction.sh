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


#### Pre-calculations for CisBP reconstructions
eukarspecies='species.dat'
selectspecies='species_handcraft_measured2.nwk.txt'

recdir=${outdir}'Reconstruction/'

# Download files from cisbp
# Main files: proteins.tab, tfs.tab, prot_features.tab, domains.tab
# create directory for every species
# run:
# python parse_cisbp.py
# fused protein sequences and names are provided in rrmkhfused.tar.gz


# All fastas for longest isoforms
domainfasta='_RRMKH_domain.fasta'
combinedfasta='_RRMKH_domain_combined.fasta'
fusedfasta='_RRMKH_domain_fused.fasta'

perform=${outdir}Performance/
#mkdir ${perform}JPLEunknown
#mkdir ${perform}LRunknown

full=0

if [ $full = 1 ]; then
cd $recdir
specieslist=$(ls -d */)
for species in $specieslist
do
# Compute sequence identity to RNAcompete experiments
python ../../Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${species}${species}${combinedfasta} --alignment local Blosum62 -11 -1 --savemax --lengthnorm alignment --outname ${species}${species}${combinedfasta%.fasta}_2Rncmmpt_lnormalignment_alignlocal_id --secondset ${outdir}Rncmpt.aaseq.ext15_combined.fasta

# Generate protein sequence features 5 posmismatch k-mers
# gapped 5-mers for JPLE
python ../../Scripts/protein-features_beta.py --mprocessing 4 --domainfasta ${species}${species}${fusedfasta} --k 5 --outname ${species}${species}${combinedfasta%.fasta} --posmismatchkmer 1
done

seqid2rncmpt='_RRMKH_domain_combined_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq.npz'
features='_RRMKH_domain_fused_5mer_features.npz'

##### Pre-calculations for Sequence search #####

#Extract protein features from 49 species for proteins that have less than 30% seq id to any protein in RNCMPT as dummy profiles
python Scripts/dummpy_proteins.py $seqid2rncmpt $features 30 $selectspecies
# Create set of dummy protein sequences for UMAP figure
python Scripts/dummy_proteins.py $seqid2rncmpt $features 95 $selectspecies
# Embed dummy protein sequences
python Scripts/Specificity_predictors.py None None --JPLE svd 122 significant_response lsq global --savelatenttest --normY2 --normP2 --recompute ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef.npz Reconstruction/Dummy_RRMKH_49species_5pm1_features_idmax95.0.npz none none JPLE_RRMKHcomplete/RecomputeDummy

fi

complete=0

if [ $complete = 1 ]; then
# Predict protein sequences for specificities that were left out in trainingset. Add dummy sequences to test out-of-distribution behaviour
for i in {0..354}
do
# Run JPLE for all leave-one out sets
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz --trainingset $i ${indir}${rrmkhlist%.list}-Zscores_420_origrncmpt_pearson0.6_testsetSingle.list --JPLE svd 122 significant_response lsq local --savetestcorrelation --savetopintersection --savetopcorrelation --savetestprofileslight --savelatentstats --determine_sequence ${recdir}Dummy_RRMKH_49species_5pm1_features_idmax30.0.npz --normP2 --normY2 --outname ${perform}JPLEunknown/jplesvd122_set${i}
# Run OLS for all leave-one out sets
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz --trainingset $i ${indir}${rrmkhlist%.list}-Zscores_420_origrncmpt_pearson0.6_testsetSingle.list --LinearRegression 0 False False --savetestcorrelation --savetopintersection --savetopcorrelation --savetestprofileslight --determine_sequence ${recdir}Dummy_RRMKH_49species_5pm1_features_idmax30.0.npz --normP2 --normY2 --outname ${perform}LRunknown/lr_set${i}
done
fi

## Interface predictions from RNAcompte experimental data
# map JPLE k-mer weights back to domain sequences
python Scripts/kmer_map2seq.py --kmerfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_PreconstructionZ.npz --original_kmerfeatures ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz --domainfasta ${outdir}Rncmpt.aaseq.ext15_domain_fused.fasta --normalize --multiprocessing 4
# Split fused domain sequences into domains again as input to visualization
python Scripts/splitdomain_importance.py ${outdir}Rncmpt.aaseq.ext15_domain.fasta ${outdir}Rncmpt.aaseq.ext15_domain_fusedjplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence.fasta > ${outdir}Rncmpt.aaseq.ext15_domain_fusedjplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta



