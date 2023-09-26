#!/usr/bin/bash -l 

function prepend_to_basename () {
    addition=$1
    filename=$2
    my_dir=$(dirname $filename)
    my_basename=$(basename $filename)
    my_new_filename="$my_dir/${addition}$my_basename"
    echo $my_new_filename
}

#echo "$( prepend_to_basename stuff_ 06-DE_lists/deseq2-BORED_vs_EXCITED.csv )"

processed_counts_file=$1
mapping_file=$2
if [ -z "$mapping_file" ]; then
    mapping_file=gene_name_mapping.txt
fi

my_bin=$(dirname $0)

tabbed_counts_file=$(prepend_to_basename tabbed $processed_counts_file)
sed 's/,/\t/g' $processed_counts_file > $tabbed_counts_file

sorted_file=$(prepend_to_basename sorted_ $tabbed_counts_file)

# Save the header
head -n1 $tabbed_counts_file > $sorted_file

# sort the rest
tail -n +2 $tabbed_counts_file | sort >> $sorted_file

# Remove tabbed file
rm $tabbed_counts_file

# Join two files using the first column
cat $sorted_file | join --header -t $'\t' -j1 $mapping_file  -
