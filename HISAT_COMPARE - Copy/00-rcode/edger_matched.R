#!/bin/env Rscript
# styler: off
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=1-00:00
# styler: on

# Load the library while suppressing verbose messages.
suppressPackageStartupMessages(library(edgeR))

# Output stream.
outfile <- ""

args <- commandArgs(trailingOnly = TRUE)

count_file <- args[1]
nameA <- args[2] # Assuming this is the Control
nameB <- args[3]

counts <- read.table(count_file, header = TRUE, sep = "\t", row.names = 1)

# Create search pattern for each group
# Assumes name then underscore
# WARNING: To be more precise, need to prepend this pattern with word boundary
patternA <- paste0(nameA, "_")
patternB <- paste0(nameB, "_")

# Get boolean list of which columns match each group
is_groupA <- grepl(patternA, colnames(counts))
is_groupB <- grepl(patternB, colnames(counts))


# Get column names for each group
samplesA <- colnames(counts[, is_groupA])
samplesB <- colnames(counts[, is_groupB])

counts <- counts[, c(samplesA, samplesB)]

# Count number of samples in each group
num_groupA <- sum(is_groupA)
num_groupB <- sum(is_groupB)

# Get list of groups corresponding to selected columns
# Force it to stay in the right order (actual name would cause alphabetical
#   sorting, which would tend to be wrong about half of the time.
condA <- rep("cond1", num_groupA)
condB <- rep("cond2", num_groupB)
group <- c(condA, condB)


# Creates a DGEList object from a table of counts and group.
dge <- DGEList(counts = counts, group = group)

# Maximizes negative binomial conditional common likelihood to estimate a
#   common dispersion value across all genes.
dis <- estimateCommonDisp(dge)

# Estimates tagwise dispersion values by an empirical Bayes method based on
#   weighted conditional maximum likelihood.
tag <- estimateTagwiseDisp(dis)

# Compute genewise exact tests for differences in the means between the groups.
etx <- exactTest(tag)

# Extracts the most differentially expressed genes.
etp <- topTags(etx, n = nrow(counts))

# Get the scale of the data
scale <- dge$samples$lib.size * dge$samples$norm.factors

# Get the normalized counts
normed <- round(t(t(counts) / scale) * mean(scale))

# Select the edger result dataframe.
data <- etp$table

# Create column placeholders.
data$baseMean <- 1
data$baseMeanA <- 1
data$baseMeanB <- 1
data$foldChange <- 2^data$logFC

# Rename the column.
names(data)[names(data) == "logFC"] <- "log2FoldChange"

# Compute the adjusted p-value
data$PAdj <- p.adjust(data$PValue, method = "hochberg")

# Reorganize the columns. Make errPct the last column before normalized data.
# ORIGINAL order
# 1 log2FoldChange
# 2 logCPM
# 3 PValue
# 4 FDR
# 5 baseMean
# 6 baseMeanA
# 7 baseMeanB
# 8 foldChange
# 9 PAdj

data <- data[c(
  "baseMean",
  "baseMeanA",
  "baseMeanB",
  "foldChange",
  "log2FoldChange",
  "logCPM",
  "PValue",
  "PAdj",
  "FDR"             # column 9
)]

# Sort the data by PValue to compute false discovery counts.
data <- data[with(data, order(PValue, -foldChange)), ]

# Compute the false discovery counts on the sorted table. # column 10
data$falsePos <- seq_len(nrow(data)) * data$FDR

# Create a merged output that contains the normalized counts. # column 11
total <- merge(data, normed, by = "row.names")

# Sort the data for the output.
total <- total[with(total, order(PValue, -foldChange)), ]

# Rename columns for consistency with other methods.
colnames(total)[1] <- "name"

# Start index for normalized data for condition 1.
start1 <- 12

# End index for normalized data for condition 1.
end1 <- start1 + num_groupA - 1

# Start index for normalized data for condition 2.
start2 <- end1 + 1

# End index for normalized data for condition 2.
end2 <- start2 + num_groupB - 1

total$baseMean <- rowMeans(total[start1:end2])
total$baseMeanA <- rowMeans(total[, start1:end1])
total$baseMeanB <- rowMeans(total[, start2:end2])

# Round the numbers
total$foldChange <- round(total$foldChange, 3)
total$FDR <- round(total$FDR, 4)
total$PAdj <- round(total$PAdj, 4)
total$logCPM <- round(total$logCPM, 1)
total$log2FoldChange <- round(total$log2FoldChange, 1)
total$baseMean <- round(total$baseMean, 1)
total$baseMeanA <- round(total$baseMeanA, 1)
total$baseMeanB <- round(total$baseMeanB, 1)
total$falsePos <- round(total$falsePos, 0)

# Write the result to the standard output.
write.csv(total, file = outfile, row.names = FALSE, quote = FALSE)
