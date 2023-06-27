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

recdir=${outdir}'Reconstruction/'
crmgdir=${outdir}'CRMGs/'

specieslist='species_handcraft_measured4.nwk.txt'
speciesnwk='species_handcraft_measured4.nwk'
# get evolutionary distance matrix from nwk
python3 Scripts/nwkdist.py ${crmgdir}${speciesnwk}
speciesdist='species_handcraft_measured4.nwkdistance_matrix.txt'

# latent representation for each species
latreps='Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent.npz'
simcut='0.2'
measlatrep=${outdir}'JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal_-testset_Platent.npz'
shortproteinid=${crmgdir}'species_handcraft_measured4_RRMKH_domains_concatenated_lnormminid.npz'
longproteinid=${crmgdir}'species_handcraft_measured4_RRMKH_domains_concatenated_lnormalignmentid.npz'
protfasta=${crmgdir}'species_handcraft_measured4_RRMKH_domains_concatenated.fasta'

# Figure 5

cd $crmgdir
# Computation of the Conserved RNA motif Groups (CRMGs)
python3 ../../Scripts/CRMG_clustering.py --analyze ${specieslist} -1 $latreps 0.2 ../../${measlatrep} ${speciesdist} ../../${shortproteinid} 70 30 --correct_domaincount --correct_length ../../${protfasta} 55 --allidlink
crmgfile=JPLEcluster_species_handcraft_measured4.nwk-1_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent_clusterscut0.2_idlink.txt
crmgclusters=JPLEcluster_species_handcraft_measured4.nwk-1_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent_clusterscut0.2_idlink_cluster.txt

# Condensed crmgfile 
python3 ../../Scripts/make_CRGfile.py ../../${crmgfile}
crmginfo=JPLEcluster_species_handcraft_measured4.nwk-1_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent_clusterscut0.2_idlink.info.txt

# File with stats for each present day species
python3 ../../Scripts/crmg_stats_perspecies.py ../../${crmgclusters}
crmgspecies=JPLEcluster_species_handcraft_measured4.nwk-1_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_Platent_clusterscut0.2_idlink_cluster_species_numbers.txt

# A and B)
# plot tree from nkw
python3 ../../Scripts/plottree.py ${speciesnwk} --outname ${fig5}${speciesnwk%nwk}.svg
# Plot CRMG evolution figures
python3 ../../CRMG_clustering.py --datasheet $crmgfile $specieslist -1 $speciesdist --savefig --plot_branchpoints --evopercentage
cd ../../

# C)
cisbpdom='Cisbp-protdomains.txt'
# Determine new measurements from genetic clustering method
python2 Scripts/new_measurements.py $latreps --removeconfident $reccosine --determine_next --domain_only $cisbpdom  RRM --savefig --logx
python2 Scripts/new_measurements.py $latreps --removeconfident $reccosine --determine_next --domain_only $cisbpdom KH --savefig --logx
# redo with provided files
python2 Scripts/new_measurements.py Recovered_proposed_next_measurement_cut0.2-domainKH.txt,Recovered_proposed_next_measurement_cut0.2-domainRRM.txt "0.2 KH,0.2 RRM" --savefig
# Find growth curve parameters to translate CRMG numbers into number of bound 7-mers by any RBP group
python3 Scripts/kmer-addition.py Data/Zscores_420_origrncmpt.txt $rpbzoolatrep 0.2 --savefig --domain_split Data/RBP_420domainclass.txt RRM,KH --rescalex 3200,93.44732647442399 439,100



# Supplementary figures

# File with cosine distance to closest measured RBP
reccosine=${outdir}${recdir}'CisBP_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat'
# file with kingdoms for all cisbp species
specking=${outdir}${recdir}'speciesclades.txt'
# Determine how many RBPs one could reconstruct with different cutoffs to a measured rbp
python2 Scripts/jple_reconstruction.py $reccosine $cisbpdom $specking --savefig

# Peptide representations of RBPs from T. vaginalis
tvagpeptides=${outdir}'Reconstruction/Trichomonas_vaginalis/Trichomonas_vaginalis_RRMKH_domain_fused_5mer_features.npz'
# JPLE model Eigenvectors
coeficients=${outdir}'JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef.npz'
# Estimate false negative rates between two unmeasured RBPs depending on similarity to measured ones
python2 Scripts/Specificity_predictors_latentvariance.py $tvagpeptides $coeficients --fraction 0.127,0.2 --meanx

# Investigate false positive rate for KH domain falling into RRMs. Cluster on JPLE only and detect fraction of clusters with both Domain classes present
python2 Scripts/new_measurements.domain.py $latreps $cisbpdom $measlatrep --determine_next --confidence 0.2
# Plot results for selected cluster thresholds
python2 Scripts/new_measurements.domain.py JPLEcosineclusters_cut0.2.txt 0.2 --savefig

# generate hard training sets from sequence identity trees generated by nearest neighbor joining 
python3 Scripts/visualize/tree_based_sets.py --seqidfile ${outdir}/Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt --domain_sensitive ${data}/RBP_420domainclass.txt --savefig ${supfig}NJ0065tree_hardset --dpi 200 --linkage nj --cutoff 0.065 --rbplist ${data}/RNCMPT_unique_experiments_RRM,KH.list
for i in {0.49}; do python Scripts/Specificity_predictors.py Zscores_420_origrncmpt.txt $trainfeatures --trainingset $i Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id_nj0.065_sets.txt  --JPLE svd 122 significant_response lsq local --savetestcorrelation --savetopintersection --savelatentstats --normP2 --normY2 --outname Performance/JPLE06/jplesvd122_set${i} ; done
python2 Scripts/calibration.py 50 JPLEnj0065/jplesvd122_set,_svdsignificant_response122.0_maplsq_declocal_-testset_latentdists.npz both cosine Outdir/Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt Data/Zscores_420_origrncmpt_pearson.mat --savefig --reconstruction_cut 0.6
# Find alternative threshold for confidence across all pairs
# latentdistances between all RBPs in embedding


# Check cluster assignment visually within neighbor joining trees for plants, vertebrates, and invertebrates for rrm and kh separately 
python3 linkage.2.py --seqidfile ${longproteinid} Micromonas_pusilla/${latreps},Physcomitrella_patens/${latreps},Selaginella_moellendorffii/${latreps},Arabidopsis_thaliana/${latreps},Cannabis_sativa/${latreps},Amborella_trichopoda/${latreps} "M.pusilla,P.patens,S.moellendorffii,A.thaliana,C.sativa,A.trichopoda" --numbermatrix --classage --dpi 60 --domain_select KH  --cluster ${crmgfile} --classage --linkage nj --outname ${supfig}species_handcraft_measured4_RRMKH_domains_concatenated_lnormalign

python3 linkage.2.py --seqidfile ${longproteinid} Micromonas_pusilla/${latreps},Physcomitrella_patens/${latreps},Selaginella_moellendorffii/${latreps},Arabidopsis_thaliana/${latreps},Cannabis_sativa/${latreps},Amborella_trichopoda/${latreps} "M.pusilla,P.patens,S.moellendorffii,A.thaliana,C.sativa,A.trichopoda" --numbermatrix --classage --dpi 60 --domain_select RRM  --cluster ${crmgfile} --classage --linkage nj --outname ${supfig}species_handcraft_measured4_RRMKH_domains_concatenated_lnormalign
python3 linkage.2.py --seqidfile $longproteinid Homo_sapiens/$latreps,Mus_musculus/$latreps,Mondelphis_domestica/$latreps,Gallus_gallus/$latreps,Xenopus_tropicalis/$latreps,Danio_rerio/$latreps "H.sapiens,M.muculus,M.domestica,G.gallus,X.tropicalis,D.rerio" --numbermatrix --classage --dpi 60 --domain_sensitive  --cluster $rssgs --classage --linkage nj --outname ${supfig}species_handcraft_measured4_RRMKH_domains_concatenated_lnormalign --domain_select RRM

python3 linkage.2.py --seqidfile $longproteinid Caenorhabtidis_elegans/$latreps,Caenorhabtidis_briggsae/$latreps,Pristionchus_pacificus/$latreps,Ascaris_suum/$latreps,Trichuris_suis/$latreps,Musca_domestica/$latreps,Drosophila_melanogaster/$latreps,Anopheles_gambiae/$latreps,Apis_mellifera/$latreps "C.elegans,C.briggsae,P.pacificus,A.suum,T.suis,M.domestica,D.melanogaster,A.gambiae,A.mellifera" --numbermatrix --classage --dpi 60 --domain_sensitive  --cluster $rssgs --classage --linkage nj --outname ${supfig}species_handcraft_measured4_RRMKH_domains_concatenated_lnormalign --domain_select RRM

python3 linkage.2.py --seqidfile $longproteinid Homo_sapiens/$latreps,Mus_musculus/$latreps,Mondelphis_domestica/$latreps,Gallus_gallus/$latreps,Xenopus_tropicalis/$latreps,Danio_rerio/$latreps "H.sapiens,M.muculus,M.domestica,G.gallus,X.tropicalis,D.rerio" --numbermatrix --classage --dpi 60 --domain_sensitive  --cluster $rssgs --classage --linkage nj --outname ${supfig}species_handcraft_measured4_RRMKH_domains_concatenated_lnormalign --domain_select KH 

python3 linkage.2.py --seqidfile $longproteinid Caenorhabtidis_elegans/$latreps,Caenorhabtidis_briggsae/$latreps,Pristionchus_pacificus/$latreps,Ascaris_suum/$latreps,Trichuris_suis/$latreps,Musca_domestica/$latreps,Drosophila_melanogaster/$latreps,Anopheles_gambiae/$latreps,Apis_mellifera/$latreps "C.elegans,C.briggsae,P.pacificus,A.suum,T.suis,M.domestica,D.melanogaster,A.gambiae,A.mellifera" --numbermatrix --classage --dpi 60 --domain_sensitive  --cluster $rssgs --classage --linkage nj --outname ${supfig}species_handcraft_measured4_RRMKH_domains_concatenated_lnormalign --domain_select KH 


cisbpcosine=${recdir}'CisBP_Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat'
# Compare the one-to-one and other orthologs with each other
python spclusters.2.py $specieslist -1 $speciesdist $latreps 0.2 0.4 $cisbpcosine 0.4 $shortproteinid --analyze --correct_length $protfasta 55

# Show evolutionary relationships with PWMs plotted in current species and connectsion between related CRMGs
python arcdiagram.py $crmgfile ../Zscores_420_origrncmpt_pwm.hmot --savefig --meanPWM -showall




