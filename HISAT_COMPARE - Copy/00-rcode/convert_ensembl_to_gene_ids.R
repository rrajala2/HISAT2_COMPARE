
setwd("C:/Users/cejdan/Dropbox (OMRF)/_Bioinformatics_training_center/Training_for_specific_labs/Lee_lab")

raw_counts = read.csv(file = "04-raw_counts/counts.txt", 
                      sep ="\t", 
                      comment.char = "#")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl")


filter_options = listFilters(ensembl)
attribute_options = listAttributes(ensembl)

ensembl_ids_table = getBM(attributes = c("ensembl_gene_id", 
                            "external_gene_name", 
                            "description"),
                          values = raw_counts$Geneid,
                          mart = ensembl)


counts_with_ids = merge(x = ensembl_ids_table, y = raw_counts, by.x = "ensembl_gene_id", by.y = "Geneid")

write.table(counts_with_ids, 
            file = "counts_geneIDs.txt", 
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)


