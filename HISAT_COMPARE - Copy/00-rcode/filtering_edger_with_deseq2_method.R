#!/bin/env Rscript

## Author: Josh Starnes (with modifications by Christopher Bottoms)
## Original source: https://statquest.org/statquest-filtering-genes-with-low-read-counts/
## NOTE: If you use this template, be sure to cite both edgeR and
## DESeq2 since this template uses DESeq2's method for identifying
## an optimal minimum sequencing depth for a gene to be
## considered transcribed.

## CITE:
## Robinson MD, McCarthy DJ and Smyth GK (2010). "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.â" Bioinformatics, 26, pp. -1.
## Love MI, Huber W and Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, pp. 550. doi: 10.1186/s13059-014-0550-8.

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(genefilter))

args = commandArgs(trailingOnly=TRUE)
count_file <- args[1] 
name_A     <- args[2] # Assuming this is the Control
name_B     <- args[3]

counts <- read.table(count_file,header=TRUE, sep="\t", row.names=1)

# Create search pattern for each group
# Assumes literal period then name then underscore
pattern_A <- paste0('\\.', name_A, '_')
pattern_B <- paste0('\\.', name_B, '_')

# Get boolean list of which columns match each group
is_group_A <- grepl(pattern_A,colnames(counts))
is_group_B <- grepl(pattern_B,colnames(counts))


# Get column names for each group
samples_A <- colnames(counts[,is_group_A])
samples_B <- colnames(counts[,is_group_B])

counts <- counts[,c(samples_A,samples_B)]

# Count number of samples in each group
num_group_A <- sum(is_group_A)
num_group_B <- sum(is_group_B)

# Get list of groups corresponding to selected columns 
# This forces it to stay in the right order. Using the actual name will cause alphabetical sorting, which may be wrong half of the time.
cond_A <- rep("cond1",num_group_A)
cond_B <- rep("cond2",num_group_B)
group <- c(cond_A, cond_B)

# Creates a DGEList object from a table of counts and group.
dge <- DGEList(counts=counts, group=group)

###
### CHANGE THE VALUE FOR "root.name".
### root.name is the the start of all output file names. For example, if you
### set root.name like this:
###
### root.name <- "wt_v_ko_edgeR_20141210_"
###
### then the output files will be named:
###
### wt_v_ko_edgeR_20141210_results.txt <- that contains all statistical tests
### wt_v_ko_edgeR_20141210_sig_results.txt <- just the genes with FDR < 0.05
### wt_v_ko_edgeR_20141210_mds.pdf <- the MDS (like PCA) plot
### wt_v_ko_edgeR_20141210_smear.pdf <- the smear plot
###
root.name <- "results_"

###
### That said, you can run the analysis without actually saving the results to
### files. If you want to do that, change the following variable
###
#save.files <- FALSE # if you don't want to actually save the output, use FALSE
save.files <- TRUE # if you want to save the output to files, use TRUE


###
### raw.data is still just a general "table" (described above). Here we are
### converting that table into something that is specific to edgeR.
y <- DGEList(
             counts=raw.data[,2:ncol(raw.data)],
             genes=raw.data[,1]
     )

###
### Here we are just doing more things that edgeR requires.
###
rownames(y) <- 1:nrow(y$counts)
y <- calcNormFactors(y)

###
### Here is where we tell edgeR how to group the samples. 1 represnts one group
### and 2 represents another group.
### CHANGE THIS IF YOU CHANGED THE FILES YOU READ IN.
### EXAMPLES:
### If you read in 3 WT and 3 KO, you want: 1,1,1,2,2,2
###
### If you read in 4 WT and 2 KO, you want: 1,1,1,1,2,2
###
### If you read in 2 WT and 6 KO, you want: 1,1,2,2,2,2,2,2
###
replicate <- factor(c(1,1,1,2,2,2))

###
### The following commands set up the tests and then execute them
###
design <- model.matrix(~replicate)
rownames(design) <- colnames(y)
design

## Notes about the next three fucntion calls...
## These can be consolodated to one call by using "estimateDisp(y, design)"
## I'm not doing that because I only want the estimated dispersion
## parameters (Disp and BCV) printed to the screen, and estimateDisp() doesn't
## give you that option. It prints a ton of stuff, or nothing at all.
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

y.fit <- glmFit(y, design)
y.lrt <- glmLRT(y.fit)

##
## Here we figure out the optimal way to filter out lowly expressed genes.
## We do this so that the FDR correction has more power - it doesn't filter
## out as many True Positives. The few genes that FDR adjusts, the more true
## positives we'll get (however, there's a point where you are filtering out
## true positives, so you can't just toss out all genes)
##
reads.cpm <- cpm(y)
unfiltered.results <- data.frame(id=y$genes, reads.cpm, y.lrt$table)

## NOTE: If you use the following code to se
## We're making "filter" equal the second largest CPM value for
## each gene, rather than the average CPM value for each gene.
## Setting it to the second largest CPM makes filter act like
## edgeR's selection where at least 2 samples have to have
## CPMs > MIN.CPM
filter <- apply(X=reads.cpm,
                MARGIN=1,   # Operate across rows
                FUN=function(data) {data[order(rank(data), decreasing=TRUE)[2]]}
          )

lowerQuantile <- mean(filter == 0) # The quantile for which the second largest CPM value is zero

# Exclude the top 5% of values, unless all non-zero values are in the top 5%.
if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1

# Theta is fifty values ranging from the lower quantile to the upper quantile
theta <- seq(lowerQuantile, upperQuantile, length=50)

filtPadj <- filtered_p(filter=filter,
                       test=unfiltered.results$PValue,
                       theta=theta,
                       method="BH")

min.fdr <- 0.05
numRej <- colSums(filtPadj < min.fdr, na.rm = TRUE)

filter.quantiles <- quantile(filter, probs=theta)

lo.fit.theta <- lowess(numRej ~ theta, f=1/5)

if (max(numRej) <= 10) {
    j <- 1
} else {
    residual <- if (all(numRej==0)) {
                    0
                } else {
                    numRej[numRej > 0] - lo.fit.theta$y[numRej > 0]
                }
    thresh <- max(lo.fit.theta$y) - sqrt(mean(residual^2))
    j <- if (any(numRej > thresh)) {
             which(numRej > thresh)[1]
         } else {
             1
         }
}

plot(theta,
     numRej,
     type="b",
     xlab="",
     ylab="",
     lwd=3,
     frame.plot=FALSE,
     col="black"
)

if (save.files) {
    file.name = paste(root.name, "cpm_thresholds.pdf", sep="")
    print(file.name)
    pdf(file=file.name)

        plot(theta,
             numRej,
             main="Threshold for optimal FDR correction",
             type="b",
             lwd=3,
             frame.plot=FALSE,
             col="black",
             ylab="# Significant Genes",
             xlab="Quantile"
        )

        lines(lo.fit.theta$x[1:21],
              lo.fit.theta$y[1:21],
              col="red",
              lwd=3
        )
        abline(h=thresh,
               col="blue",
               lwd=3
        )

    dev.off()
}

filtered.results <- unfiltered.results
filtered.results$FDR <- filtPadj[, j, drop=TRUE]

filtered.results$de <- sign(filtered.results$logFC)*(filtered.results$FDR < 0.05)
filtered.results$sig <- abs(filtered.results$de)

filtered.results[is.na(filtered.results$sig),]$sig <- 0

###
### This command does a "head" on just the significantly DE genes. This let's
### us, at a glance, make sure things worked as expected.
###
head(filtered.results[filtered.results$sig==1,])

draw_final_plot <- function(insig_data,sig_data) {
    plot(x=insig_data$logCPM,
         y=insig_data$logFC,
         pch=16,
         col="#00000011",
         ylab="logFC",
         xlab="logCPM",
         ylim=c(y.min, y.max),
         xlim=c(x.min, x.max)
    )
    points(x=sig_data$logCPM,
           y=sig_data$logFC,
           pch=20,
           col="#FF000088"
    )
}

if (save.files) {
    ###
    ### save results
    ###
    file.name = paste(root.name, "results.txt", sep="")
    print(file.name)
    write.table(filtered.results, file=file.name, sep="\t", row.names=FALSE, quote=FALSE)
    
    file.name = paste(root.name, "sig_results.txt", sep="")
    print(file.name)
    write.table(filtered.results[filtered.results$sig==1,], file=file.name, sep="\t", row.names=FALSE, quote=FALSE)
    
    ###
    ### save mds plot
    ###
    file.name = paste(root.name, "mds.pdf", sep="")
    print(file.name)
    pdf(file=file.name)
        plotMDS(y)
    dev.off()
    
    ###
    ### save smear plot
    ###
    file.name = paste(root.name, "smear.pdf", sep="")
    print(file.name)
    pdf(file=file.name, width=8, height=4)
        y.min <- min(filtered.results$logFC)
        y.max <- max(filtered.results$logFC)
        if (abs(y.min) > y.max) {
            y.max = abs(y.min)
        } else {
            y.min = -1 * y.max
        }
        x.min <- min(filtered.results$logCPM)
        x.max <- max(filtered.results$logCPM)
        insig <- subset(filtered.results, FDR > 0.05)
        sig <- subset(filtered.results, FDR <= 0.05)
        plot(x=insig$logCPM,
             y=insig$logFC,
             pch=16,
             col="#00000011",
             ylab="logFC",
             xlab="logCPM",
             ylim=c(y.min, y.max),
             xlim=c(x.min, x.max)
        )
        points(x=sig$logCPM,
               y=sig$logFC,
               pch=20,
               col="#FF000088"
        )
    dev.off()
}

# Plot same plot in interactive viewer
plot(x=insig$logCPM,
     y=insig$logFC,
     pch=16,
     col="#00000011",
     ylab="logFC",
     xlab="logCPM",
     ylim=c(y.min, y.max),
     xlim=c(x.min, x.max)
)
points(x=sig$logCPM,
       y=sig$logFC,
       pch=20,
       col="#FF000088"
)
