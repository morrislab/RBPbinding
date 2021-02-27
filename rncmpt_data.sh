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


# similarity of setA and setB
reproducible='Fullset_z_scores-setA-vs-setB-comparison.txt' # somehow has only 422? # something was wrong with renormalizing and lost e-scores from original normalization

# Folder for all slide shows
slides='Slideshows'


### General data manipulation for Figure1 and subsequent analsysis 

# calculate motif identity matrices
python Scripts/Motif_identity_matrices.py ${indir}${zscores} --zscores --outname ${outdir}${zscores%.txt}_pearson.mat --pearson --savetxt
python Scripts/Motif_identity_matrices.py ${indir}${zscores} --zscores --outname ${outdir}${zscores%.txt}_top10.mat --topnum 10 --topset --savetxt
python Scripts/Motif_identity_matrices.py ${indir}${zscores} --zscores --outname ${outdir}${zscores%.txt}_top100.mat --topnum 100 --topset --savetxt

# Generate PWMs for all measured specificities
python Scripts/extract_motifs_zscores.py --zscores ${indir}${zscores} --outname ${outdir}${zscores%.txt}_pwm highpwm --kmerchoice top 10 --clusterparams 1 3 1 4
# Extract IUPACs
python ../motifsextract.py Zscores_420_origrncmpt_pwm.hmot > Zscores_420_origrncmpt_pwm.motif
# plot PWMs
python Scripts/pylogo.py ${outdir}${zscores%.txt}_pwm.hmot

# extract domain fastas
python Scripts/make_extdomain_fasta.py ${indir}${masterfile} 0 ${outdir}Rncmpt.aaseq

# extract domain fastas with 15 AA flanking sequences
python Scripts/make_extdomain_fasta.py ${indir}${masterfile} 15 ${outdir}Rncmpt.aaseq.ext15


# compute sequece identity on long protein sequence
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${outdir}Rncmpt.aaseq.ext15_combined.fasta --alignment local Blosum62 -11 -1 --savetxt --lengthnorm alignment --outname ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id
# compute protein sequence identities on short sequence
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${outdir}Rncmpt.aaseq.ext15_combined.fasta --alignment local Blosum62 -11 -1 --savetxt --lengthnorm min --outname ${outdir}Rncmpt.aaseq.ext15_combined_lnormmin_alignlocal_id
# compute sequence identities for highest domain
python /Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${outdir}Rncmpt.aaseq_domain.fasta --alignment local Blosum62 -11 -1 --savetxt --lengthnorm alignment --outname ${outdir}Rncmpt.aaseq_domain_lnormalignment_alignlocal_id
python Scripts/Fuse_domain_identity_from_multseq_aligment.py --domainidentities ${outdir}Rncmpt.aaseq_domain_lnormalignment_alignlocal_id.txt --ident_option highdom

# clustering by protein sequence
# long sequence norm
python3 Scripts/agglomerative_clustering.py ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt sim 30 single gt
# short sequence norm
python3 Scripts/agglomerative_clustering.py ${outdir}Rncmpt.aaseq.ext15_combined_lnormmin_alignlocal_id.txt sim 40 single gt
# clustering by specificity
python3 Scripts/agglomerative_clustering.py ${outdir}${zscores%.txt}_pearson.mat sim 0.6 complete gt
# clustering sequence clusters by specificity
python3 Scripts/agglomerative_clustering.py ${outdir}${zscores%.txt}_pearson.mat sim 0.6 complete gt --clusters ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id_aggclusters-singlegt30.npz


