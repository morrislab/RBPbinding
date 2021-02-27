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




#### Pre-computations for Figure 4
# Reconstruct specificities for species protein sequences with JPLE
cd Reconstruction
specieslist=$(ls -d */)
for species in $specieslist
do
python Scripts/Specificity_predictors.py ${indir}${zscores} ${outdir}Rncmpt.aaseq.ext15_domain_fused_5mer_features.npz --JPLE svd 122 significant_response lsq local --saveclosestcorrelation --savelatentstats --savetestprofileslight --savetestlatent --normY2 --normP2 --recompute ${outdir}JPLE_RRMKHcomplete/jplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef.npz ${species}${species%/}_RRMKH_domain_fused_5mer_features.npz local none ${species}Reconstruction
done


# Reconstruct protein specificities from sequence identity
for species in $specieslist
do
# Reconstruct protein specificities with specificities from  Nature 2013 paper only
python Scripts/topid.py ${species}${species%/}_RRMKH_domain_combined_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq.npz --maxexperiment 291
# Nature 2013 + this study
python Scripts/topid.py ${species}${species%/}_RRMKH_domain_combined_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq.npz
done

# Generate PWMs for reconstructed specificities
for species in $specieslist
do
cd $l
python Scripts/extract_motifs_zscores.py --predictlight Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_profiles_light.txt --outname Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal_pwm highpwm --kmerchoice top 10 --clusterparams 1 3 1 4 --normalize_zscores
cd ..
done

# Get CisBP statistics
# JPLE cut-off 0.127, sequence identity cutoff 70% and 40%
# Nature 2013 + this study, with precise cutoff seqid
python Scripts/cisbp-reconst-stats.2.py Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-directfit_highcorrelation.dat  _RRMKH_domain_combined_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq_maxid.txt _proteinnames.txt 0.127 70 ../Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6.npz cisbp2
# Nature 2013 + this study, with tolerant cutoff seqid
python Scripts/cisbp-reconst-stats.2.py Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-directfit_highcorrelation.dat  _RRMKH_domain_combined_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq_maxid.txt _proteinnames.txt 0.127 40 ../Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6.npz cisbp2
# Nature 2013, with precise cutoff seqid
python Scripts/cisbp-reconst-stats.2.py Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-testset_latentstats.dat Reconstructionjplesvd122_svdsignificant_response122.0_maplsq_decglobal__coef_recdeclocal-directfit_highcorrelation.dat  _RRMKH_domain_combined_2Rncmmpt_lnormalignment_alignlocal_id.on.Rncmpt.aaseq-maxexp291_maxid.txt _proteinnames.txt 0.127 70 ../Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6.npz cisbp06 --maxmeasured 291
cisbpreconstructionstats='CisBP_reconstruction-stats-jple0.127-id70.0-cisbp2-Zscores_420_origrncmpt_pearson_aggclusters-completegt0.6_nspecies691.csv'

