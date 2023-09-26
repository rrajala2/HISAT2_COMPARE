#!/bin/awk

# Nested ternary operator (see https://en.wikipedia.org/wiki/%3F: )
function notempty_or(x,y,z) {
    return x?x:y?y:z
}

# Split lines by tabs
BEGIN{FS="\t"}

{
    # remove semicolon at the end of line
    sub(/;$/,"")

    # "globally" remove every double-quote characters from attributes column
    gsub(/"/,"",$9)

    # Clear variables from previous iterations 
    for (i in attr) { delete attr[i] }
    gene_id=""
    gene_name=""

    # Skip lines that aren't a "gene" or "exon" feature
    if($3!~"gene" && $3!~"exon") { next }

    # Split the attribute column into attribute strings 
    split($9, a, /;\s*/) # (splitting on semicolon with optional space)

    # Split each attribute string into its key value pair
    for (i in a) {
        split(a[i], j, / |=/) # Space delimiter for GTF, equal sign for GFF3
        attr[ j[1] ] = j[2]   # Store key-value attribute pair
    }

    # Extract gene ID
    gene_id = notempty_or(attr["gene_id"],attr["ID"],"NA")

    # Skip records without a gene_id
    if(gene_id ~ /^NA$/) {next}

    # Remove "gene:" prefix if gene_id came from "ID" field
    if(gene_id == attr["ID"]) {
        sub(/gene:/, "",gene_id)
    }

    # Skip transcript records
    if(gene_id ~ /transcript/) {next}

    # Get gene name (or default to gene_id)
    gene_name = notempty_or(attr["gene_name"],attr["Name"],gene_id)

    # Print tab-delimited gene ID and name
    print gene_id "\t" gene_name 
}
