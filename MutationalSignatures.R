# Load Packages 
library(BSgenome) # for available ref-genome 
library(MutationalPatterns) # for mutational signature analysis 
library(NMF) # for non-negative matrix factorization(NMF or NNMF)
library(gridExtra) # for managing plots and grid system 

# List vailable genomes 
head(available.genomes())

# Load data 
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Store VCF Files into vcf_files[vcf=variant call format] 
vcf_files <- list.files(system.file("extdata", package = "MutationalPatterns"), pattern = ".vcf", full.names = TRUE)
vcf_files

# Sample names [came from vcf_files]
sample_names <- c("colon1", "colon2", "colon3", "intestine1", "intestine2", "intestine3", 
                  "liver1", "liver2", "liver3")

# Load VCF files into GRanges[vcf_files, sample_names, ref_genome]
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

# 96 Mutational Profile 
mut_mat <- mut_matrix(vcf_list=vcfs, ref_genome = ref_genome)
head(mut_mat)


# De novo mutational signature extraction using NMF
mut_mat <- mut_mat + 0.0001 
estimate <- nmf(mut_mat, rank=2:5, method = "brunet", nrun=10, seed=123456) 
plot(estimate) 

# Heatmap 
consensusmap(estimate)

# Mutational Signatures Extraction 
nmf_res <- extract_signatures(mut_mat, rank = 2, nrun=10)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B")
plot_96_profile(nmf_res$signatures, condensed = TRUE)

# Contribution of Signatures in a bar plot 
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signatures, mode = "relative")
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signatures, mode = "absolute")
grid.arrange(pc1, pc2)

# Flip X and Y coordinates 
pc3 <- plot_contribution(nmf_res$contribution, nmf_res$signatures, mode="relative", coord_flip = TRUE)
pc4 <- plot_contribution(nmf_res$contribution, nmf_res$signatures, mode="absolute", coord_flip = TRUE)
grid.arrange(pc3, pc4)

# Heatmaps 
plot_contribution_heatmap(nmf_res$contribution, sig_order = c("Signature A", "Signature B"))
plot_contribution_heatmap(nmf_res$contribution, cluster_samples = FALSE)

# GRID 
pch1 <- plot_contribution_heatmap(nmf_res$contribution, sig_order = c("Signature A", "Signature B"))
pch2 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples = FALSE)
grid.arrange(pch1, pch2, widths=c(2, 1.6))

# Compare reconstructed mutational profile with original mutational profile 
plot_compare_profiles(mut_mat[,1], 
                      nmf_res$reconstructed[,1], 
                      profile_names = c("Original", "Reconstructed"), 
                      condensed = TRUE)


# Cosmic Mutational Signatures 
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", sep = " ")

cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standardnew_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)

cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

# Plot mutational profile of the first two COSMIC signatures
plot_96_profile(cancer_signatures[,1:2], condensed = TRUE, ymax = 0.3)

# Dendogram 
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
plot(hclust_cosmic)

cos_sim(mut_mat[,1], cancer_signatures[,1])
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
plot_cosine_heatmap(cos_sim_samples_signatures,
                    col_order = cosmic_order,
                    cluster_rows = TRUE)

