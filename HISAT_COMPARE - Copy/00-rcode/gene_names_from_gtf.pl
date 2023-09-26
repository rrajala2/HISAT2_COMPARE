#!/bin/env perl
use strict;
use warnings;

use feature 'say';

my %gene_name_for;

while( my $line = readline) {

    # Remove newline
    chomp $line;

    # Skip comments (lines that start with hashtag)
    next if $line =~ /\A \# /xms; 

    # replace double quotes with nothing
    $line =~ s/"//gxms;
    
    # Attributes are all in the ninth field
    my $attributes = (split /\t/, $line)[8];

    # break up attribute pairs by semicolon delimiter 
    my @attribute_strings = (split /\s*;\s*/, $attributes);

    my %val_for;
    for my $attrib_string (@attribute_strings) {

        # Split attribute strings using one or more spaces as delimiter
        # WARNING: If there is a space in the gene name, we're in trouble
        my ($key, $val) = split /\s+/, $attrib_string;

        # Skip if key is empty
        next if ! $key;

        # Skip if val is empty
        die "Missing value for $key" if ! $val;

        # Store value
        $val_for{$key} = $val;
    }

    # Store info if gene id is found
    if (exists $val_for{"gene_id"}) {
        my ($gene_id, $gene_name) = @val_for{"gene_id", "gene_name"};

        # For empty gene names, just use the gene id
        $gene_name ||= $gene_id;

        # Warn about overwriting different gene name
        if (exists $gene_name_for{$gene_id} &&
                $gene_name_for{$gene_id} ne $gene_name) {
            warn "Previous gene_name for $gene_id is $gene_name_for{$gene_id}. Replacing it with $gene_name";
        }

        # Store gene name for gene id
        $gene_name_for{$gene_id} = $gene_name; 
    }
}

# Get an ordered list of gene ids
my @ids = sort(keys %gene_name_for);

# Print out the gene id to name mapping
for my $id (@ids) {
    say "$id\t$gene_name_for{$id}";
}
