#!/usr/bin/bash -l
# Remove comment lines | # Remove bam directory and extension from sample names
grep -v -P '^#' $@     | sed 's/\t03-bam\/\(\w\+\).bam/\t\1/g'
