#!/usr/bin/bash -l 

# Get the directory that this script is in
bin_dir=$(dirname $0)

my_input=$1

my_cat=cat

# If file is compressed, then use zcat instead of cat
if [[ "$my_input" =~ gz ]]; then
    my_cat=zcat
fi

echo -e "name\tgene_name"

$my_cat $my_input                            | `# stream file`                        \
    awk -f $bin_dir/gene_names_from_gtf.awk  | `# extract gene names with awk script` \
    sort -u                                    `# sort and only keep unique entries`

# backslash at the end of a line is a "line continuation", which means the current line
#   and next line are effectively one long line
# inline quotes are of the format `# comments `
#   (Everything inside backticks is evaluated and the result interpolated. Since comments
#    are ignored, this returns nothing and has no effect. But the comments are still
#    useful to us humans.)
