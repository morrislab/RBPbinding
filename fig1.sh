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


seqid=${outdir}'Rncmpt.aaseq.ext15_combined_lnormalignment_alignlocal_id.txt'
mpearson=${zscores%.txt}'_pearson.mat'

# piechart figure 1
python Scripts/piechart.py ${indir}${masterfile} ${fig1}

### Figure 1A # Motif similarity matrix after sequence identity clustering
python Scripts/make_newick_phylo.motifs.3.py $seqid single None 'Seq identity' 100.,50.,0. ${indir}${masterfile} ${outdir}${zscores%.txt}_pwm.hmot --savefig $fig1 150 --showsimilarity $mpearson --markset 291 --assignspecies --proteinlist ${indir}${highestzscoreset}

### Figure 1B # Motif similarity matrix after sequence identity clustering for a-rich cluster
python Scripts/make_newick_phylo.motifs.3.py $seqid single None 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig $fig1 150 --showsimilarity $mpearson --markset 291 --assignspecies --proteinlist ${indir}pearsonlist-arich.txt

# Figure 1C # Individual z-score scatter plots for sart3 family
python Scripts/make_newick_phylo.motifs.3.py $seqid single 30 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig ${fig1} 200 --showindividual ${indir}${zscores} --proteinlist ${indir}pearsonlist-sart3.list

# Figure 1D Sequence identity to Pearson correlation relationship
python Scripts/pearson_vs_id.py 1 $mpearson $seqid --protlist ${indir}${rrmkhlist} --boxplot --ylabel 'Pearson' --selfcorr Fullset_z_scores-setA-vs-setB-comparison.txt 1 --setylim 0 1




####### Supplementary figures 1 #####
## Figure SX. Motif similarity matrix after sequence identity clustering for chosen clusters
lists='pearsonlist1.txt pearsonlist2.txt pearsonlist3.txt pearsonlist4.txt pearsonlist5.txt pearsonlist6.txt pearsonlist7.txt pearsonlistgacrich.txt pearsonlistgrich.txt pearsonlisturich.txt'
for l in $lists
do
python Scripts/make_newick_phylo.motifs.3.py $seqid single None 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig $fig1 150 --showsimilarity $mpearson --markset 291 --assignspecies --proteinlist ${indir}${l}
done

## Figure SX all individual z-score scatter plots
python Scripts/make_newick_phylo.motifs.3.py ${outdir}Rncmpt.aaseq_domain_lnormalignment_alignlocal_id_combined_highdom.txt single 30 'Seq identity' 100.,50.,0. $masterfile ${outdir}${zscores%.txt}_pwm.hmot --savefig ${fig1} 200 --showindividual ${indir}${zscores}

#Figure SX Motif similarity clustering
python Scripts/make_newick_phylo.motifs.3.py $mpearson complete 0.6 'Pearson' 1.,0.5,0. ${indir}${masterfile} ${outdir}${zscores%.txt}_pwm.hmot --savefig '' 100 --plot --markset 291 --assignspecies --proteinlist ${indir}${highestzscoreset}
python Scripts/make_newick_phylo.motifs.3.py $mpearson complete None 'Pearson' 1.,0.5,0. ${indir}${masterfile} ${outdir}${zscores%.txt}_pwm.hmot --savefig '' 100 --plot --markset 291 --assignspecies --proteinlist ${indir}${highestzscoreset}

# Sequence identity to motif similarity for KH, RRM and top100
python Scripts/pearson_vs_id.py 1 $mpearson $seqid --protlist ${indir}RNCMPT_unique_experiments_KH.list --boxplot --ylabel 'Pearson' --selfcorr ${indir}${reproducible} 1 --setylim 0 1
python Scripts/pearson_vs_id.py 1 $mpearson $seqid --protlist ${indir}RNCMPT_unique_experiments_RRM.list --boxplot --ylabel 'Pearson' --selfcorr ${indir}${reproducible} 1 --setylim 0 1
python Scripts/pearson_vs_id.py 1 ${outdir}${zscores%.txt}_top100.mat $seqid --protlist ${indir}RNCMPT_unique_experiments_RRM,KH.list --boxplot --ylabel 'Pearson' --selfcorr ${indir}${reproducible} 3 --setylim 0 1


