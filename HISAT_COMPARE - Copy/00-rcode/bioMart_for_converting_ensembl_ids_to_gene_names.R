library(biomaRt)

default_id_column_name <- "name"

# Command line argument.
args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
    id_column_name = args[1]
} else {
    id_column_name = default_id_column_name
}

de_list = read.csv(file = "stdin",
                      comment.char = "#")

ensembl <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl") #NOTE: To find other datasets to use:
                                                         #          library(biomaRt)
                                                         #          listDatasets(ensembl)

# List options

# Restrict output based on start and end position, chrom name, gene name, description...
# filter_options = listFilters(ensembl)

# attributes are what you want to pull out of the database
# attribute_options = listAttributes(ensembl)

ensembl_ids_table = getBM(attributes = c("ensembl_gene_id",
                                         "external_gene_name",
                                         "description"),
                          values = de_list[[id_column_name]],
                          mart = ensembl)


counts_with_ids = merge(x = ensembl_ids_table, y = de_list, by.x = "ensembl_gene_id", by.y = id_column_name)

# Write to stdout
write.table(counts_with_ids,
            file = "",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
