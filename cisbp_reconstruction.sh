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

reconst='Outdir/Reconstruction/'

# Generate domain fasta containing only measured RBPs with RRMs and KH domains
python Scripts/extractRRMKH.py ${indir}$rrmkhlist ${outdir}Rncmpt.aaseq.ext15_domain.fasta > ${outdir}Rncmpt.aaseq.ext15_domain_RRMKH.fasta
full=0

#### Generate files for CisBP #### 
# Run interface reconstructions
list=$(ls -d $reconst)
# Alternatively, only create content for list of spcies
#list=("Outdir/Reconstruction/Arabidopsis_thaliana/")

if [ full = 1 ];then
for l in $list; do echo $l; python Scripts/Specificity_predictors.py none none --JPLE svd 122 significant_response lsq global --save_reconstructionPtoP --normY2 --normP2 --recompute ${outdir}JPLE_RRMKHcomplete/jplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef.npz ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_fused_5mer_features.npz global none ${l}Reconstruction ; done

# Map peptide weights onto domain sequences for confident reconstructions (jple dist < 0.127) to protein sequence 
for l in $list
do
python Scripts/kmer_map2seq.py --kmerfeatures ${l}Reconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP.npz --original_kmerfeatures ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_fused_5mer_features.npz --domainfasta ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_fused.fasta --normalize --multiprocessing 4
# split interface masks into individual domains
python Scripts/splitdomain_importance.py ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain.fasta ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_fusedReconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP-norm_kmer2sequence.fasta > ${l}Reconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP-norm_kmer2sequence_domain.fasta
# Compute domain identities to RNCMPT measurements
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain.fasta --alignment local Blosum62 -11 -1 --savemax --lengthnorm alignment --outname ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_2Rncmmpt_lnormalignment_alignlocal_id --secondset ${outdir}Rncmpt.aaseq.ext15_domain_RRMKH.fasta
# Give domains secondary structure profiles of measured proteins
python Scripts/align_secondarystructure_template.py ${outdir}RRM_Imask/Rncmpt.RRMKH.ss.split.rn.fasta ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain.fasta ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq-max.txt
done

complete=0
for l in $list
do
# Make jpgs with interface scores and location of domain in sequences for confident reconstructions (jple dist < 0.127) to protein sequence
mkdir ${l}Interfacescores
if [ $complete = 1 ]; then
python Scripts/interface_secstruc_map.py ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_secondary_templates.fasta Reconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP-norm_kmer2sequence_domain.fasta --fullseqloc ${l}"$(echo $l | cut -d'/' -f 3 )"_proteinsequences_longestisoform.fasta --ylimmean --savefig ${l}Interfacescores/JPLEscore.jpg --confidence Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat 0.127
else
python Scripts/interface_secstruc_map.py ${l}"$(echo $l | cut -d'/' -f 3 )"_RRMKH_domain_secondary_templates.fasta Reconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP-norm_kmer2sequence_domain.fasta --ylimmean --savefig ${l}Interfacescores/JPLEscore.jpg --confidence Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat 0.127
fi
done

for l in $list
do
# Generate PWMs
python Scripts/extract_motifs_zscores.py --predictlight ${l}Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_profiles_light.txt --outname ${l}Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal_pwm highpwm --kmerchoice top 10 --clusterparams 1 3 1 4 --normalize_zscores

# Filter for confident top100 profiles
python Scripts/confidentzscore.py ${l}Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat ${l}Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_profiles_light.txt 0.127
done
fi

for l in $list
do
mkdir ${l}PWMs
# Filter for confident PWM profiles
python Scripts/confidentpwm.py ${l}Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat ${l}Reconstructionjplesvd122_svdsignificant_response122.hmot 0.127
python Scripts/pylogo.py ${l}Reconstructionjplesvd122_svdsignificant_response122_confidence0.127.hmot --infocont --removeframe --outdir ${l}PWMs/
done











