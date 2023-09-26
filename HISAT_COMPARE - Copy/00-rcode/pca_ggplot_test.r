
# Load the library while suppressing verbose messages.


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))


# Output stream (empty because it will be to standard out instead of an actual file).
outfile = ""

# Set the plot dimensions.
WIDTH = 12
HEIGHT = 13

args = commandArgs(trailingOnly=TRUE)

count_file <- args[1] # Counts file name needs to be supplied as command line argument

# Just for testing locally, delete this line:
# count_file = "//data/cbds/hpc-nobackup/bitbucket/rna-seq/level_1_rnaseq_pipeline/04-raw_counts/counts.txt"

counts <- read.table(count_file,header=TRUE, sep="\t", row.names=1)


matrix_counts = 
counts %>%
  select(!contains(c("Geneid", "Chr", "Start", "End", "Strand", "Length"))) %>% # removes the matched columns
  mutate(across(where(is.integer), as.numeric)) %>% # changes all integers to numeric for pca
  rename_with( ~ str_replace(.x,
                             ".*\\.bam\\.(.*)\\.bam",
                             "\\1"))
  #filter( ~ sum(.x) != 0)

matrix_counts = matrix_counts[rowSums(matrix_counts[])>0,]
  
test_scale = as.data.frame(scale(matrix_counts, center = TRUE, scale = TRUE))



pca = prcomp(matrix_counts, center = TRUE, scale = TRUE) # Need to scale for PCA to cluster accurately

# Commented out plot (which defaults to creating file named "Rplots.pdf")
# plot(pca$rotation)


pca_df = as.data.frame(pca$rotation)

max_pc1 = max(pca_df$PC1)
min_pc1 = min(pca_df$PC2)

max_pc2 = max(pca_df$PC2)
min_pc2 = min(pca_df$PC2)

axis_max = max(0.3, max_pc1, max_pc2) + 0.2
axis_min = min(-0.3, min_pc1, min_pc2) - 0.2

pc1_percent_var = round(pca[["sdev"]][1] / sum(pca[["sdev"]]),4)*100
pc2_percent_var = round(pca[["sdev"]][2] / sum(pca[["sdev"]]),4)*100
#
#library(devtools)
##install_github("vqv/ggbiplot")
#library(ggbiplot)
#
#ggbiplot(pca,circle=T,obs.scale=1,varname.size=3)
#ggbiplot(pca,circle=T,choices=c(1,3),obs.scale=1,varname.size=20)
#

# Open the drawing device (i.e. standard out as a PDF file).
pdf(file ='|cat', width = WIDTH, height = HEIGHT)

    ggplot(data = pca_df, aes(x = PC1, y = PC2, label = row.names(pca_df))) +
      geom_point() +
      geom_text_repel() +
      coord_fixed() +
      scale_x_continuous(name=paste0("PC1 (", pc1_percent_var, "% of var)"), limits=c(axis_min, axis_max)) +
      scale_y_continuous(name=paste0("PC2 (", pc2_percent_var, "% of var)"), limits=c(axis_min, axis_max)) +
      ggtitle("Principal Components 1 and 2")
    
    
    
#    # Draw the PCA plots using base R. Could improve the look with ggplot or other packages.
#    plot(pca_df$PC2 ~ pca_df$PC1, xlab = "PC1", ylab = "PC2", main = "PCA") +
#    text(pca_df$PC2 ~ pca_df$PC1, labels = row.names(pca_df),cex=0.5, font=1)

# Turn off the device.
dev.off()
