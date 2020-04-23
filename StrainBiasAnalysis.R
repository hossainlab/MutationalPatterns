# Load Packages 
library(BSgenome) # for available ref-genome 
library(MutationalPatterns) # for mutational signature analysis 
library(NMF) # for non-negative matrix factorization(NMF or NNMF)
library(gridExtra) # for managing plots and grid system 
library("TxDb.Hsapiens.UCSC.hg19.knownGene")


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


# Set tissue type [came from vcf_files]
tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))



# Load Genes 
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes_hg19

# Mutation Strand 
strand = mut_strand(vcfs[[1]], genes_hg19)
head(strand, 10)

# Transcriptional strand bias analysis
# Mutation Strand Matrix 
mut_mat_strand <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
# Matrix Indexing 
mut_mat_strand[1:5, 1:5]

# Starin Counts 
strand_counts <- strand_occurrences(mut_mat_strand, by=tissue)
head(strand_counts)

# Strain Bias 
strand_bias <- strand_bias_test(strand_counts)
strand_bias

# Plotting 
ps1 <- plot_strand(strand_counts, mode = "relative")
ps2 <- plot_strand_bias(strand_bias)
grid.arrange(ps1, ps2)


# Replicative strand bias analysis
repli_file = system.file("extdata/ReplicationDirectionRegions.bed",
                         package = "MutationalPatterns")

# Read Replication Direction Regions File
repli_strand = read.table(repli_file, header = TRUE)
# Store in GRanges Object

repli_strand_granges = GRanges(seqnames = repli_strand$Chr,
                            ranges = IRanges(start = repli_strand$Start + 1,
                            end = repli_strand$Stop),
                            strand_info = repli_strand$Class)

# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) = "UCSC"
repli_strand_granges

strand_rep <- mut_strand(vcfs[[1]], repli_strand_granges, mode = "replication")
head(strand_rep, 10)

# Create Mutation Strand Matrix 
mut_mat_strad_rep <- mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges,
                                       mode = "replication")
# Matrix Indexing 
mut_mat_strain_rep[1:5, 1:5]

repli_strand_granges$strand_info <- factor(repli_strand_granges$strand_info,
                                           
                                             levels = c("right", "left"))

mut_mat_s_rep2 <- mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges,
                                  
                                        mode = "replication")
mut_mat_s_rep2[1:5, 1:5]

strand_counts_rep <- strand_occurrences(mut_mat_strain_rep, by=tissue)
head(strand_counts)
strand_bias_rep <- strand_bias_test(strand_counts_rep)
strand_bias_rep


ps1 <- plot_strand(strand_counts_rep, mode = "relative")
ps2 <- plot_strand_bias(strand_bias_rep)
grid.arrange(ps1, ps2)


# Extract signatures with strand bias
nmf_res_strand <- extract_signatures(mut_mat_strand, rank = 2)
colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")
a <- plot_192_profile(nmf_res_strand$signatures, condensed = TRUE)
b <- plot_signature_strand_bias(nmf_res_strand$signatures)
grid.arrange(a, b, ncol = 2, widths = c(5, 1.8))
