# Load Packages 
library(BSgenome) # for available ref-genome 
library(MutationalPatterns) # for mutational signature analysis 
library(NMF) # for non-negative matrix factorization(NMF or NNMF)
library(gridExtra) # for managing plots and grid system 
library(biomaRt)

chromosomes <- seqnames(get(ref_genome))[1:22]
plot_rainfall(vcfs[[1]], title = names(vcfs[1]),
              chromosomes = chromosomes, cex = 1.5, ylim = 1e+09)


regulatory <- useEnsembl(biomart="regulation",
             dataset="hsapiens_regulatory_feature", GRCh = 37)

CTCF <- getBM(attributes = c('chromosome_name',
                            'chromosome_start',
                            'chromosome_end',
                            'feature_type_name',
                            'cell_type_name'),
                             filters = "regulatory_feature_type_name",
                            values = "CTCF Binding Site",
                              mart = regulatory)
