#!/usr/bin/env Rscript
#
# Draws a heatmap from the output that contains a normalized matrix
#
library(stringr)

# Set the plot dimensions.
WIDTH = 12
HEIGHT = 13

# Command line argument.
args = commandArgs(trailingOnly=TRUE)

# Optional argument for FDR cutoff.
if (length(args) > 0) {
  LIMIT = as.integer(args[1])
} else {
  # Default FDR cutoff
  LIMIT = 5
}

# Input values should be in percent!
LIMIT = LIMIT/100

# Set the margins
MARGINS = c(9, 12)

# Relative heights of the rows in the plot.
LHEI = c(1, 5)

# Load the library.
suppressPackageStartupMessages(library(gplots))

# Read normalized counts from the standard input.
data = read.csv("stdin", header=TRUE, as.is=TRUE, sep="\t")

count_under_limit <- sum(data$FDR <= LIMIT, na.rm=TRUE)

if (count_under_limit < 1) {
    print(paste0("No data is below the FDR rate of ",LIMIT))

    # Quit because there is no data to plot
    quit(status=1)
}

# Subset data for values under a treshold.
data = subset(data, data$FDR <= LIMIT)

# The heatmap row names will be extracted from the first columns.
row_names = paste0(
                data[, 1], # Ensemble ID
                "(",

                str_trim(
                    paste0(
                        data[,2],  # more common name
                        " ",
                        data[,3]  # description
                    )
                ),

                ")"
            )

# The normalized data starts past the rightmost of these columns.
idx = which(colnames(data) == "falsePos") + 1

# The normalized counts are on the right size.
counts = data[, idx : ncol(data)]

# Load the data from the second column on.
values = as.matrix(counts)

# Adds a little noise to each element to avoid the
# clustering function failure on zero variance rows.
values = jitter(values, factor = 1, amount = 0.00001)

# Normalize each row to a z-score
zscores = NULL
for (i in 1 : nrow(values)) {
    row = values[i,]
    zrow = (row - mean(row)) / sd(row)
    zscores = rbind(zscores, zrow)
}

# Set the row names on the zscores.
row.names(zscores) = row_names

# Turn the data into a matrix for heatmap2.
zscores = as.matrix(zscores)

# Open the drawing device.
pdf('|cat', width = WIDTH, height = HEIGHT)

    # Set the color palette.
    color_palette = greenred

    # Clean up column names
    sample_names <- colnames(zscores)

    sample_names <- gsub('bam\\.|\\.bam|X03\\.','',sample_names)

    # Draw the heatmap.
    heatmap.2(
        zscores,
        col=color_palette,
        labCol=sample_names,
        density.info="none",
        Colv=NULL,
        dendrogram="row",
        trace="none",
        margins=MARGINS,
        lhei=LHEI
    )

# Finalize and close the PDF file
dev.off()

