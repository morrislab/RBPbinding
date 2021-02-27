

#### Generate files for CisBP #### 
# Run interface reconstructions
for l in $list; do echo $l;python Scripts/Specificity_predictors.py none none --JPLE svd 122 significant_response lsq global --save_reconstructionPtoP --normY2 --normP2 --recompute ${outdir}JPLE_RRMKHcomplete/jplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef.npz ${l}${l%/}_RRMKH_domain_fused_5mer_features.npz global none ${l}Reconstruction ; done

# Map peptide weights onto domain sequences for confident reconstructions (jple dist < 0.127) to protein sequence 
for l in $list
do
cd $l
python Scripts/kmer_map2seq.py --kmerfeatures Reconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP.npz --original_kmerfeatures ${l%/}_RRMKH_domain_fused_5mer_features.npz --domainfasta ${l%/}_RRMKH_domain_fused.fasta --normalize --multiprocessing 4

# split interface masks into individual domains
python Scripts/splitdomain_importance.py ${l%/}_RRMKH_domain.fasta ${l%/}_RRMKH_domain_fusedReconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP-norm_kmer2sequence.fasta > Reconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP-norm_kmer2sequence_domain.fasta
# Compute domain identities to RNCMPT measurements
for l in $list; do python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${l}${l%/}_RRMKH_domain.fasta --alignment local Blosum62 -11 -1 --savemax --lengthnorm alignment --outname ${l}${l%/}_RRMKH_domain_2Rncmmpt_lnormalignment_alignlocal_id --secondset ../Rncmpt.aaseq.ext15_domain_RRMKH.fasta; done
# Give domains secondary structure profiles of measured proteins
python Specificity/align_secondarystructure_template.py ${indir}Rncmpt.RRMKH.ss ${l%/}_RRMKH_domain.fasta ${l%/}_RRMKH_domain_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq-max.txt
# Make jpgs with interface scores and location of domain in sequences for confident reconstructions (jple dist < 0.127) to protein sequence
mkdir Interfacescores
python Scripts/interface_secstruc_map.py ${l%/}_RRMKH_domain_secondary_templates.fasta Reconstructionjplesvd0.95hom_svdsignificant_response0.95_maplsq_decglobal__coef_recdecglobal-testset_PreconstructionP-norm_kmer2sequence_domain.fasta --fullseqloc ${l%/}_proteinsequences_longestisoform.fasta --ylimmean --savefig Interfacescores/JPLEscore.jpg --confidence Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat 0.127

# Generate PWMs
python Scripts/extract_motifs_zscores.py --predictlight Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_profiles_light.txt --outname Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal_pwm highpwm --kmerchoice top 10 --clusterparams 1 3 1 4 --normalize_zscores
python Scripts/confidentpwm.py ${l}Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat ${l}Reconstructionjplesvd122_svdsignificant_response122.hmot 0.127
python Scripts/confidentzscore.py Arabidopsis_thaliana/Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat Arabidopsis_thaliana/Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_profiles_light.txt 0.127

cd ..
done











