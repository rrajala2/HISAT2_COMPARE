
# Load the library while suppressing verbose messages.

if (!require(dplyr)) {
  install.packages("dplyr", dependencies = TRUE, quiet = TRUE)
}

suppressPackageStartupMessages(library(dplyr))


# Output stream (empty means standard out, instead of an actual file).
outfile <- ""

# Set the plot dimensions.
plot_width <- 12
plot_height <- 13

args <- commandArgs(trailingOnly = TRUE)

count_file <- args[1] # first command line argument

counts <- read.table(count_file, header = TRUE, sep = "\t", row.names = 1, as.is=TRUE)


matrix_counts <-
  counts %>%
  select(!contains(c(
    "Geneid", "Chr", "Start", "End", "Strand", "Length", # featureCount columns
    "target_id", "length", "eff_length", "X1", "X2", "X3" # kallisto-related columns
  ))) %>% # removes the matched columns
  mutate(across(where(is.integer), as.numeric)) # changes all integers to numeric for pca


pca <- prcomp(matrix_counts, scale = TRUE) # Need to scale for PCA to cluster accurately

# Commented out plot (which defaults to creating file named "Rplots.pdf")
# plot(pca$rotation)


pca_df <- as.data.frame(pca$rotation)

# Open the drawing device (i.e. standard out as a PDF file).
pdf(file = "|cat", width = plot_width, height = plot_height)

# Draw the PCA plots using base R. Could improve look with ggplot or other pkg.
plot(pca_df$PC2 ~ pca_df$PC1, xlab = "PC1", ylab = "PC2", main = "PCA") +
  text(pca_df$PC2 ~ pca_df$PC1, labels = row.names(pca_df), cex = 0.5, font = 1)

# Turn off the device.
dev.off()
