library(readr)
output_ggCmanSeqCluster_FindGenes <- read_csv("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC/2024_03_12/output_ggCmanSeqCluster_FindGenes.txt")

View(output_ggCmanSeqCluster_FindGenes)
colnames(output_ggCmanSeqCluster_FindGenes) <- output_ggCmanSeqCluster_FindGenes[10,1]
output_ggCmanSeqCluster_FindGenes <- output_ggCmanSeqCluster_FindGenes[-(1:10),]

likelihoods_FindGenes_df <- data.frame(matrix(0, ncol = 4, nrow = 50))
#output_ggCmanSeqCluster_FindGenes[9,1]
colnames(likelihoods_FindGenes_df) <- c("ggCmanSeqCl","COGmanSeqCl","ggCPP","COGPP")

for (i in 1:49) {
  likelihoods_FindGenes_df[i,1] <- strsplit(output_ggCmanSeqCluster_FindGenes[9 * i,1]$`[1] "ggCaller_manSeqClusters"`, "\\[1\\]")[[1]][2]
}

#sort.list(likelihoods_FindGenes_df$ggCmanSeqCl) # should be decreasing=TRUE but it does not seem to recognise the minus sign
sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]


# Find the genes that are in these categories:
intermed_gene_presence_absence_consensus <- readRDS(file = "ggC_intermed_gene_presence_absence_consensus.rds")
# compute boolean gene vectors
n_groups <- ceiling((nrow(intermed_gene_presence_absence_consensus)-1)/25)
find_genes_df <- data.frame(matrix(0,nrow = nrow(intermed_gene_presence_absence_consensus)-1, ncol = 50))
for (i in 1:25) {
  find_genes_df[,i] <- rep(c(rep(0,24),1),(n_groups + 4))[i:(nrow(intermed_gene_presence_absence_consensus)-2 + i)]
  find_genes_df[,i+25] <- c(rep(0,(nrow(intermed_gene_presence_absence_consensus) - 1-n_groups)),rep(1,n_groups),rep(0,(nrow(intermed_gene_presence_absence_consensus) -n_groups)))[(1 + (i-1) * n_groups) : (nrow(intermed_gene_presence_absence_consensus)-1 + (i-1) * n_groups)]
}

rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])==1) # 571
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])==2) # 69

new_fitting_vec <- rep(0, nrow(intermed_gene_presence_absence_consensus)-1)
new_fitting_vec <- as.integer(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])>=2)
saveRDS(new_fitting_vec, "ggCmanSeqCl_FindGenesResults1.rds")

# Fit for ggCmanSeq 30 44 46 21  5 49 47 22  1 43 (all genes that appear at least twice)

#sigma_f           m           v 
#0.176375633 0.009003682 0.061537212 
#[1] "det_mcmc_2 log likelihood"
#log_likelihood 
#-260.6447 
#[1] "det_mcmc_2 mean log likelihood"
#[1] -260.1266
### That is a relatively good, but not outstanding, fit for the 3-param model.
# there was one fit with a log-likelihood of -253.9187 (likelihoods_FindGenes_df$ggCmanSeqCl[30])

new_fitting_vec2 <- as.integer(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])>=1)
saveRDS(new_fitting_vec2, "ggCmanSeqCl_FindGenesResults2.rds")

#    sigma_f          m          v 
#0.06377337 0.02299571 0.18618946 
#[1] "det_mcmc_2 log likelihood"
#log_likelihood 
#-261.6746 
#[1] "det_mcmc_2 mean log likelihood"
#[1] -258.4692

### ggCPP FindGenes version
output_ggCPP_FindGenes <- read_csv("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC/2024_03_12/output_ggCPP_FindGenes.txt")
View(output_ggCPP_FindGenes)
colnames(output_ggCPP_FindGenes) <- output_ggCPP_FindGenes[10,1]
output_ggCPP_FindGenes <- output_ggCPP_FindGenes[-(1:10),]
for (i in 1:49) {
  likelihoods_FindGenes_df[i,3] <- strsplit(output_ggCPP_FindGenes[9 * i,1]$`[1] "ggCaller_PopPUNK"`, "\\[1\\]")[[1]][2]
}
#sort.list(likelihoods_FindGenes_df$ggCmanSeqCl) # should be decreasing=TRUE but it does not seem to recognise the minus sign
sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]
#  5 44 46 21 12 38 30  2 19 28

rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]])
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]])==1) # 568
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]])==2) # 71


# 23.07.2025
ga_ppxsero <- readRDS("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_07_23/GeneticAlg/gann.rds")
ga_ppxsero <- readRDS("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_07_28/GeneticAlg/gann.rds")
plot(ga_ppxsero)
plot(ga_ppxsero, ylim = c(-310, -200))
abline(h = -214, col = "red", lty = "dashed")
abline(h = -215, col = "red", lty = "dashed")
abline(h = -261, col = "red", lty = "dashed")
abline(h = -271, col = "red", lty = "dashed")

delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
ga_ppxsero_NFDS_vec <- as.vector(t(apply(ga_ppxsero@solution, 1, decode2)))
names(ga_ppxsero_NFDS_vec) <- names(delta_ranking)
sum(ga_ppxsero_NFDS_vec)/1934 # 0.4172699

ga_ppxsero_NFDS_vec_delta_sorted <- ga_ppxsero_NFDS_vec[names(sort(delta_ranking))]
plot(ga_ppxsero_NFDS_vec_delta_sorted)
mean(ga_ppxsero_NFDS_vec_delta_sorted[1:floor(1934 * 0.35142922)])
mean(ga_ppxsero_NFDS_vec_delta_sorted[-(1:floor(1934 * 0.35142922))])
plot(cumsum(ga_ppxsero_NFDS_vec_delta_sorted))
abline(a=0,b=mean(ga_ppxsero_NFDS_vec))
abline(v = floor(1934 * 0.35142922), col = "red", lty = "dashed")

plot(cumsum(ga_ppxsero_NFDS_vec_delta_sorted) - mean(ga_ppxsero_NFDS_vec) * (1:1934))
abline(v = floor(1934 * 0.35142922), col = "red", lty = "dashed")

# check function of those genes
ggCaller_us_rowname_genename_dict <- readRDS("ggCaller_us_rowname_genename_dict.rds")
bakta_annotation_panaroo <- read.delim("~/Documents/PhD_Project/Data/StrepPneumo_UKUSNepal/bakta_annotation/create_tsv_annotation/bakta_annotation_and_before.tsv")
bakta_annotation_panaroo_dict <- bakta_annotation_panaroo[,1]
names(bakta_annotation_panaroo_dict) <- bakta_annotation_panaroo[,2]
# functions of genes, under NFDS according to genetic algorithm:
bakta_annotation_panaroo_dict[ggCaller_us_rowname_genename_dict[names(which(ga_ppxsero_NFDS_vec ==1))]]
unname(bakta_annotation_panaroo_dict[ggCaller_us_rowname_genename_dict[names(which(ga_ppxsero_NFDS_vec ==1))]])

# and under NFDS according to prop_f
bakta_annotation_panaroo_dict[ggCaller_us_rowname_genename_dict[names(which(delta_ranking <= 0.35142922 * length(delta_ranking)))]]
unname(bakta_annotation_panaroo_dict[ggCaller_us_rowname_genename_dict[names(which(delta_ranking <= 0.35142922 * length(delta_ranking)))]])

delta_underNFDS <- rep(0, length(delta_ranking))
names(delta_underNFDS) <- names(delta_ranking)
delta_underNFDS[names(which(delta_ranking <= 0.35142922 * length(delta_ranking)))] <- 1


# overlap delta and genetic algorithm
# expected number of genes: 0.4172699 * 0.35142922 * 1934 = 283
overlap_delta_genAlg <- intersect(names(which(delta_underNFDS ==1)), names(which(ga_ppxsero_NFDS_vec==1)))
length(overlap_delta_genAlg)

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(names(which(delta_underNFDS ==1)), names(which(ga_ppxsero_NFDS_vec==1))),
  category.names = c("delta", "genAlg"),
  filename = '../venn_diagramm_NFDSgenes_delta_genAlg.png',
  output=TRUE, 
  
  # Output features
  imagetype="png" ,
  height = 960 , 
  width = 960 , 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#B3E2CD", "#FDCDAC"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

# 30.07.2025
ga_ppxsero_Nepal <- readRDS("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_07_30/GeneticAlg_NepalUK/FindGenes_Nepal_gann.rds")
ga_ppxsero_UK <- readRDS("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_07_30/GeneticAlg_NepalUK/FindGenes_UK_gann.rds")
plot(ga_ppxsero_Nepal)
abline(h = -1130.253, col = "red", lty = "dashed") # 4-param

plot(ga_ppxsero_UK, ylim = c(-700, -450))
abline(h =  -463.2552469, col = "red", lty = "dashed") # 4-param

