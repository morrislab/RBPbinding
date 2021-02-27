masterfile='S1.In_Progress.v6.For_all.DR-2.csv'
# transform Zscore npz into txt file to read and read in 
python transformnpz.py
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



##### Pre-calculations for interface study #####
## Prepare pdbs, generate interface masks
pdbss='ssrbp.txt' ## file with secondary structure profiles for rrm and kh proteins (download ss.txt.gz from pdb)
rrmpdbs=$(ls ${indir}pdbs/RRM/*pdb.gz)

mkdir ${outdir}RRM_Imask # directory for interface mask fastas
Imask='${outdir}RRM_Imask'

cd $Imask
# compute all interface masks
# extract interface residues, bound RNA motifs, the type of interaction, backbone, sugar, residue, multiple-interactions, and distance to nucleotides
for r in $rrmpdbs; do python Scripts/get_inteface_beta.py $r --interdistcut 5; done


## combine fasta files for downstream manipulation
python Scripts/combineint.py lt5.int RRMpdb
cat *fasta > RRMpdb.fasta
# search for Pfam domains in all pdb sequences
hmmsearch -o RRMpdb.rrm1.out ${indir}RRM_1.hmm RRMpdb.fasta
# extract domain fasta with 35 AA offset to combine all RRM pairs with their linkers included (different HMM from what was used for CisBP, ~20 AA shorter)
python Scripts/Hmmer_extract_domainfasta.py RRMpdb.rrm1.out RRMpdb.fasta 0.1 35 RRM
pdbrrmdomain='RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.fasta'
# concatenate domains that directly touch
python ../../concatenatedomain.py RRMpdb_Eval0.1-ext35_RRMpdb.rrm1.info
pdbrrmconcat='RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated.fasta'
# compute identity between pdb structures to filter for unique structures
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences $pdbrrmconcat --outname RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated.id --alignment local Blosum62 -11 -1 --savetxt --lengthnorm min

## Find masks for RNCMPT experiment
# Compute identity to extended RRM domains in RNCMPT
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences $pdbrrmconcat --outname RRMpdb2rncmpt.fasta.id --alignment local Blosum62 -11 -1 --savetxt --lengthnorm l1 --secondset ${outdir}Rncmpt.aaseq.ext15_domain.fasta
# Find best experiment for unique pdbs. 
# pdb masks are combined if more than 70% identical
# Experiment with more than 50% seq id and more than 3.5 motif match are chosen, and experiment with highest Id and highest motif match are used for evaluation
python Scripts/find_templates.py RRMpdb2rncmpt.fasta.id.on.Rncmpt.aaseq.npz 50 RRMpdb_motifslt5.int ${indir}Zscores_420_origrncmpt_pwm.motif 3.5 $rrmlist RRMpdblt5.int ${outdir}Rncmpt.aaseq.ext15_domain.fasta $pdbss RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated.fasta.npz $proteinnames Rncmpt_templates
# final interface masks for RNAcompete experiments
dertemplates='Rncmpt_templates-id50-msim3.5.txt'
cd ..

## Interface predictions from RNAcompte experimental data
# map JPLE k-mer weights back to domain sequences
python Scripts/kmer_map2seq.py --kmerfeatures ${outdir}JPLE_RRMKHcomplete/jplesvd0.95hom_svdsignificant_response0.95_maplsq_-testset_PreconstructionZ.npz --original_kmerfeatures ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz --domainfasta ${indir}Rncmpt.aaseq.ext15_domain_fused.fasta --normalize --multiprocessing 4
jpleint='Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_-testset_PreconstructionZ-norm_kmer2sequence.fasta'

## Visualize predictions with JPLE
# We ran scratch to predict secondary structure profiles for measured RBPs in RNAcompete
#~/SCRATCH-1D_1.2/bin/run_SCRATCH-1D_predictors.sh  Rncmpt.aaseq.ext15_domain.fasta Rncmpt.aaseq.ext15_domain.fasta.scratch 4
# File with predicted secondary structure profiles
rncmptsecstruc='Rncmpt.RRMKH.ss.split.rn.fasta'

# Split fused domain sequences into domains again as input to visualization
python Scripts/splitdomain_importance.py ${indir}Rncmpt.aaseq.ext15_domain.fasta $jpleint > Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95_nonorm_svdsignificant_response0.95_maplsq_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta
jpleintsplit='Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta'
# Generate jpg with importance map on secondary structure as CisBP resource
python Scripts/interface_secstruc_map.py ${rncmptsecstruc} ${jpleintsplit} --fullseqloc ${indir}Rncmpt.aaseq_pseq.fasta --ylimmean --savefig Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_PrecZ.jpg

#Generate pdbs with jple score for RNAcompete measures with template
python Scripts/colorpdb.py ${Imask}/4ed5_B_lt5-interface.pdb,B,C Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta,RNCMPT00121 0,1 .2 --normscore --meanlowscore --maskpdb
python Scripts/colorpdb.py ${Imask}/4c4w_E_lt5-interface.pdb,E,A Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta,RNCMPT00071 0 .2 --normscore --meanlowscore --maskpdb
python Scripts/colorpdb.py ${Imask}/2m8d_B_lt5-interface.pdb,B,A Rncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.fasta,RNCMPT00163 1 .2 --normscore --meanlowscore --maskpdb
# Show pdb in pymol with coloring by JPLE score
pymol RNCMPT00121_to_4ed5_B_CRncmpt.aaseq.ext15_domain_fusedjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal_-testset_PreconstructionZ-norm_kmer2sequence_split.pdb ${Imask}showpdb.pml
# showpdb.pml colors interface accordingly, not fully automated


# Precalculations to compare JPLE score to other measures (Can be skipped and only final files be used)
# i.e Conservation, RF with PSSM, conservation, physico-chemico features in window 5
# RRM sequences with 15AA extension from CisBP
cisbprrm='CisBP_RRM.ext15_domain_concatenated.fasta'
# search cisbp for domain homologs for pdb sequences
phmmer -o RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated.cisbp.homologs.out -E 1e-6 --domE 1e-6 --incE 1e-6 $pdbrrmconcat $cisbprrm
# extract sequences of homologs
python ~/RBPs/Scripts/extract_data/Hmmerout_getorthologs.py RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated.cisbp.homologs.out $pdbrrmconcat 0.95 0.3 $pdbrrmconcat
# Calulate conservation scores
# create individual multiple sequence alignment file
mkdir ${Imask}TemplateMSA
cd ${Imask}TemplateMSA
python Scripts/multiseq.py $pdbrrmconcat RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain_homologs_maxsim0.95_minhom0.3.fasta
cd ..
# calculate amino acid conservation with: Capra JA and Singh M. Predicting functionally important residues from sequence conservation. Bioinformatics
cd conservation_code
./conservation.sh ../TemplateMSA/  _RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain_homologs_maxsim0.95_minhom0.3.fasta
cd ..
# extract and combine conservatioin
cd ${Imask}TemplateMSA
python Scripts/extractconservation.py _RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain_homologs_maxsim0.95_minhom0.3.cons
cd ..
pdbconservation=${Imask}'TemplateMSA/Conservation_RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated_homologs_maxsim0.95_minhom0.3.txt'

# Generate ML features for pdb sequences with pssm, conservation, physicochemico features and onehot encoding
# make pssm for ML model features
cd $Imask
Scripts/pssm.py $pdbrrmconcat RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-concatenated_homologs_maxsim0.95_minhom0.3.fasta --seqweighting cutoff 0.6 --countgap --infocontent
# Generate ML features from sequence, pssm and conservation
python Scripts/extract_features.py RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.fasta --onehot --physicochemico --pssm RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain_homologs_maxsim0.95_minhom0.3cutoff0.6-gaps-ictrans-pssm.npz --conservation $pdbconservation
cd ..
pdbmlfeat='RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain-onehot-physchem-PSSMRRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain_homologs_maxsim0.95_minhom0.3cutoff0.6-gaps-ictrans-pssm-conservation-features.npz'

# Generate interface files for individual domain to train ML model on
python Scripts/domaininterface.py RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.fasta RRMpdb.int
pdbdomainint='RMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.int'

# Caculate identities between template domains to exlude "same" domain from training sets (30% domain id cut off)
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.fasta --outname RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.fasta.id --alignment local Blosum62 -11 -1 --savetxt --lengthnorm alignment --savemax
# generate leave-one-out training sets for ML model that 
python Scripts/reconstruction_trainingsets.py RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.list None None RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.fasta.id.txt 30 False Single
# Generate random forest interface predictions by leave-one-out cross-validation
# Training proteins possess less than 30% seqid to test protein
# random forest classifier, on sum of features in window of size 5
mkdir ${Imask}RF
cd ${Imask}RF
for i in {0..154}; do python Scripts/rbpinterface.py --mprocessing 4 --trainmodel $pdbdomainint --kmerfile $pdbmlfeat --trainingset ${i} RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain-RRMpdb_Eval0.1-ext35_RRMpdb.rrm1-domain.fasta.id30.0_testsetSingle.list --classifier RF 5 sum None None --outname RF/rf_set${i}.txt; done
cd ..

# Map interface predictor scores onto the RNAcompete pdb masks to get final interfaces score files
python Scripts/alignscoretomask.py $jpleint 2 $dertemplates
python Scripts/alignscoretomask.py $pdbconservation 2 $dertemplates --domainfasta
python Scripts/alignscoretomask.py RF/rf_complete_wz5sum-Nonexnorm-testinterface.fasta 3 $dertemplates --domainfasta



