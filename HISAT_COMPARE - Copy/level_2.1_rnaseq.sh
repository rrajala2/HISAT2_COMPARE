#!/usr/bin/bash -l
# Original code by István Albert
# Modifications by Christopher Bottoms and/or Nic Cejda
# Disclaimer: This code is for instructional purposes only.
# Run this thus:
#   sbatch level_2_rnaseq.sh
#SBATCH --cpus-per-task=12  # For this "level", try one cpu core per sample
#SBATCH --mem-per-cpu=6G  # For this "level", try 6G of mem per cpu core
#SBATCH --time=1-00:00
#SBATCH --mail-type END,FAIL

# Load the required modules:
module load hisat2       # for aligning to genome, output to .sam file
module load samtools     # for processing .sam files (convert to .bam and index .bam)
module load subread      # for featureCounts
module load R            # for DE-analysis
module load bioconductor # for DE-analysis
module load fastqc       # for Quality Control
module load multiqc      # for Quality Control

# Input the file path to your Indexed Genome:
GENOME="00-genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GENOME_FEATURE_FILE="00-genomes/Mus_musculus.GRCm39.107.gtf"

# Make some extra directories to help file organization:
mkdir 01-qc
mkdir 02-sam
mkdir 03-bam
mkdir 04-raw_counts
mkdir 05-pca
mkdir 06-DE_lists
mkdir 07-heatmaps

# Get gene-name mapping
sbatch --wrap="00-rcode/gene_names_from_gtf.sh $GENOME_FEATURE_FILE > gene_name_mapping.txt"

# Step 1 - Quality Control

# fastqc runs an analysis for each read.fq.gz file
fastqc --threads=12 00-reads/*.fastq.gz --outdir 01-qc

# multiqc aggregates the analysis from all the fastqc outputs.
multiqc 01-qc --outdir 01-qc

# Index the reference genome. Needs to be done only once and can be reused.
sbatch --mem=24G --cpus-per-task=12 --wait --wrap="hisat2-build -p 12 $GENOME $GENOME"


# Step 2 - align reads to genome.
cat ids.txt | parallel "hisat2 -x $GENOME -1 00-reads/{}_R1.fastq.gz   -2 00-reads/{}_R2.fastq.gz   -S 02-sam/{}.sam"

# Step 3 - sort your aligned reads. Output is a .bam file.
cat ids.txt | parallel 'samtools sort 02-sam/{}.sam > 03-bam/{}.bam'

# Step 4 - use your genome feature file to count how many reads mapped to each gene.
featureCounts -p -a $GENOME_FEATURE_FILE -M --fraction -o 04-raw_counts/counts.txt 03-bam/*.bam


# Step 5 - Check if your samples cluster together with Principal Component Analysis (PCA):

Rscript 00-rcode/pca.r 04-raw_counts/counts.txt > 05-pca/pca.pdf


# Step 6 - Run differential expression analysis with R.
#          Here are the example scripts used by the Biostar handbook for doing DE analysis.
#          Their experimental design is 3 control samples compared to 3 experiment samples. They are
#          testing if there are DE genes between the experiment group and control group.

# METHOD 1 for DE - deseq2

cat comparisons.csv | parallel --colsep=, --skip-first 'Rscript 00-rcode/deseq2.R 04-raw_counts/counts.txt {1} {2} > 06-DE_lists/deseq2-{1}_vs_{2}.csv'

# METHOD 2 for DE - edgeR use like so:

cat comparisons.csv | parallel --colsep=, --skip-first 'Rscript 00-rcode/edger.R 04-raw_counts/counts.txt {1} {2} > 06-DE_lists/edger-{1}_vs_{2}.csv'

# We can also generate some heatmaps based on the DE gene lists:

for base in $(cat comparisons.csv | parallel --colsep=, --skip-first 'echo "{1}_vs_{2}"' ); do
  cat 06-DE_lists/deseq2-$base.csv | Rscript 00-rcode/heatmap_from_csv.r > 07-heatmaps/heatmap_${base}_deseq2.pdf
  cat 06-DE_lists/edger-$base.csv | Rscript 00-rcode/heatmap_from_csv.r > 07-heatmaps/heatmap_${base}_edger.pdf

  00-rcode/add_gene_names.sh 06-DE_lists/deseq2-${base}.csv > 06-DE_lists/deseq2-${base}_w_gene_names.tsv
  00-rcode/add_gene_names.sh 06-DE_lists/edger-${base}.csv > 06-DE_lists/edger-${base}_w_gene_names.tsv
  
  cat 06-DE_lists/deseq2-${base}_w_gene_names.tsv | Rscript 00-rcode/HEATMAP_new.R > 07-heatmaps/heatmap_deseq2_${base}_w_gene_names_w_gene_names.pdf
  cat 06-DE_lists/edger-${base}_w_gene_names.tsv | Rscript 00-rcode/HEATMAP_new.R > 07-heatmaps/heatmap_edger_${base}_w_gene_names_w_gene_names.pdf
  


done  
# Get Ensembl ID to gene name mapping file









