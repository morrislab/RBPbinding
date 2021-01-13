
masterfile='S1.In_Progress.v6.For_all.DR-2.csv'
zscores='Zscores_420_origrncmpt.npz' ## --> transform into .txt file

rrmkhlist='RNCMPT_unique_experiments_RRM,KH.list'
khlist='RNCMPT_unique_experiments_KH.list'
rrmlist='RNCMPT_unique_experiments_RRM.list'

highestzscoreset='RNCMPT_unique_experiments_zcoresmax381.list'

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
# Kate needs to collect raw data, normalized data, and combine, and put in folder and compute Pearson between setA and B
reproducible='Fullset_z_scores-setA-vs-setB-comparison.txt' # somehow has only 422? # something was wrong with renormalizing and lost e-scores from original normalization






### General data manipulation for Figures

# calculate motif identity matrices
python Scripts/Motif_identity_matrices.py ${indir}${zscores} --zscores --outname ${outdir}${zscores%.txt}_pearson.mat --pearson --savetxt
python Scripts/Motif_identity_matrices.py ${indir}${zscores} --zscores --outname ${outdir}${zscores%.txt}_top10.mat --topnum 10 --topset --savetxt
python Scripts/Motif_identity_matrices.py ${indir}${zscores} --zscores --outname ${outdir}${zscores%.txt}_top100.mat --topnum 100 --topset --savetxt

# Generate PWMs for all measured specificities
python Scripts/extract_motifs_zscores.py --zscores ${indir}${zscores} --outname ${outdir}${zscores%.txt}_pwm highpwm --kmerchoice top 10 --clusterparams 1 3 1 4
# plot PWMs
python Scripts/pylogo.py ${outdir}${zscores%.txt}_pwm.hmot

# extract domain fastas
python Scripts/make_extdomain_fasta.py ${indir}${masterfile} 0 ${outdir}Rncmpt.aaseq

# extract domain fastas with 15 AA flanking sequences
python Scripts/make_extdomain_fasta.py ${indir}${masterfile} 15 ${outdir}Rncmpt.aaseq.ext15


# compute long protein sequence identities
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${outdir}Rncmpt.aaseq.ext15_combined.fasta --alignment local Blosum62 -11 -1 --savetxt --lengthnorm alignment --outname ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id
# compute short protein sequence identities
python Scripts/pairwise_alignment.2.py --mprocessing 4 --sequences ${outdir}Rncmpt.aaseq.ext15_combined.fasta --alignment local Blosum62 -11 -1 --savetxt --lengthnorm min --outname ${outdir}Rncmpt.aaseq.ext15_combined_lnormmin_alignlocal_id
# compute single domain sequence identities
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




##### Figure 1 ###### 
# piechart figure 1
python piechart.py $[indir]${masterfile} ${fig1}

### Figure 1A # Motif similarity matrix after sequence identity clustering
python Scripts/make_newick_phylo.motifs.3.py ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt single None 'Seq identity' 100.,50.,0. ${indir}${masterfile} ${outdir}${zscores%.txt}_pwm.hmot --savefig $fig1 150 --showsimilarity ${zscores%.txt}_pearson.mat --markset 291 --assignspecies --proteinlist ${indir}${highestzscoreset}

### Figure 1B # Motif similarity matrix after sequence identity clustering for a-rich cluster
python Scripts/make_newick_phylo.motifs.3.py ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt single None 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig $fig1 150 --showsimilarity ${zscores%.txt}_pearson.mat --markset 291 --assignspecies --proteinlist ${indir}pearsonlist-arich.txt

# Figure 1C # Individual z-score scatter plots for sart3 family
python Scripts/make_newick_phylo.motifs.3.py ${outdir}Rncmpt.aaseq_domain_lnormalignment_alignlocal_id_combined_highdom.txt single 30 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig ${fig1} 200 --showindividual ${indir}${zscores} --proteinlist ${indir}pearsonlist-sart3.list

# Figure 1D Sequence identity to Pearson correlation relationship
python Scripts/pearson_vs_id.py 1 ${outdir}${zscores%.txt}_pearson.mat Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt --protlist RNCMPT_unique_experiments_RRM,KH.list --boxplot --ylabel 'Pearson' --selfcorr Fullset_z_scores-setA-vs-setB-comparison.txt 1 --setylim 0 1




####### Supplementary figures 1 #####
## Figure SX. Motif similarity matrix after sequence identity clustering for chosen clusters
lists='pearsonlist1.txt pearsonlist2.txt pearsonlist3.txt pearsonlist4.txt pearsonlist5.txt pearsonlist6.txt pearsonlist7.txt pearsonlistgacrich.txt pearsonlistgrich.txt pearsonlisturich.txt'
for l in $lists
do
python Scripts/make_newick_phylo.motifs.3.py ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt single None 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig $fig1 150 --showsimilarity ${zscores%.txt}_pearson.mat --markset 291 --assignspecies --proteinlist ${indir}pearsonlist-arich.txt
done

## Figure SX all individual z-score scatter plots
python Scripts/make_newick_phylo.motifs.3.py ${outdir}Rncmpt.aaseq_domain_lnormalignment_alignlocal_id_combined_highdom.txt single 30 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig ${fig1} 200 --showindividual ${indir}${zscores}

#Figure SX Motif similarity clustering
python Scripts/make_newick_phylo.motifs.3.py ${outdir}${zscores%.txt}_pearson.mat complete 0.6 'Pearson' 1.,0.5,0. ${indir}${masterfile} ${outdir}${zscores%.txt}_pwm.hmot --savefig '' 100 --plot --markset 291 --assignspecies --proteinlist ${indir}${highestzscoreset}
python Scripts/make_newick_phylo.motifs.3.py ${outdir}${zscores%.txt}_pearson.mat complete None 'Pearson' 1.,0.5,0. ${indir}${masterfile} ${outdir}${zscores%.txt}_pwm.hmot --savefig '' 100 --plot --markset 291 --assignspecies --proteinlist ${indir}${highestzscoreset}

# Sequence identity to motif similarity for KH, RRM and top100
python Scripts/pearson_vs_id.py 1 ${outdir}${zscores%.txt}_pearson.mat ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt --protlist ${indir}RNCMPT_unique_experiments_KH.list --boxplot --ylabel 'Pearson' --selfcorr ${indir}${reproducible} 1 --setylim 0 1
python Scripts/pearson_vs_id.py 1 ${outdir}${zscores%.txt}_pearson.mat ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt --protlist ${indir}RNCMPT_unique_experiments_RRM.list --boxplot --ylabel 'Pearson' --selfcorr ${indir}${reproducible} 1 --setylim 0 1
python Scripts/pearson_vs_id.py 1 ${outdir}${zscores%.txt}_top100.mat ${outdir}Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt --protlist ${indir}RNCMPT_unique_experiments_RRM,KH.list --boxplot --ylabel 'Pearson' --selfcorr ${indir}${reproducible} 3 --setylim 0 1











