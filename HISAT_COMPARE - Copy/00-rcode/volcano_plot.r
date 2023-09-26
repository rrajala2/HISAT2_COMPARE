
#####################
# Load in Libraries #
#####################
library(EnhancedVolcano)


##############################
# Set Command Line Arguments #
##############################
outfile = ""
args = commandArgs(trailingOnly=TRUE)

file_name       <- args[1] # DE Results file
name_control    <- args[2] # Assuming this is the Control name
name_experiment <- args[3] # Experiment group name


####################################
# Read in Data and set parameters  #
####################################
de_results = read.csv(file_name)


num_top_genes = 10
FCcutoff = 2.0
pCutoff = 0.05
control = name_control
case = name_experiment
top_genes = head(de_results[order(de_results$PAdj),], num_top_genes)$name



#####################
# EnhancedVolcano   #
#####################


pdf(file = '|cat', width = 8.5, height = 11) #Set for portrait mode, Canvas is sized at (8.5 inch x 11 inch)

EnhancedVolcano(de_results,
                lab = de_results$name,
                x = 'log2FoldChange',
                y = 'PAdj',
                selectLab = top_genes,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P Adj'),
                pCutoff = pCutoff,
                FCcutoff = FCcutoff,
                pointSize = 3.0,
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabels = c("NS", 
                                 bquote(Log[2] ~ "FC " > .(FCcutoff)), 
                                 bquote("PAdj" < .(pCutoff)), 
                                 bquote("PAdj" < .(pCutoff) ~ "and" ~ Log[2] ~ "FC" > .(FCcutoff))),
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                lengthConnectors = 3,
                title= paste0(control, " vs ", case),
                subtitle = NULL
                )

dev.off()