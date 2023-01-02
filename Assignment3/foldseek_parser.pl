#!/usr/bin/perl

##Assignment 3: Perl
##Name: Kajol Tanesh Shah
##A20496724

use strict;
use warnings;

## open an output file 
open OUT, ">", 'hits_per_protein.txt';

##takes multiple file names from the @ARGV and iterate through each .fseek file one after the other 
for my $i (0..$#ARGV) {
    ##opens (reads from) the current .fseek file
    open IN, "<", $ARGV[$i];
    ## initializing a counter;
    my $count = 0; 

    ##iterate through file's content line by line
    while (my $line = <IN>) {
        chomp $line;        
        ##function split(), and grab the following columns: query, target and e-value
        my @col = split(/\t/, $line);
        ##check if the E-value is smaller or equal to 1e-20
        if ($col[10] <= 1e-20){
            ##prints the query, target and e-value (separated by semi-colons) to the hits_per_protein.txt output file
            print OUT "$col[0]; $col[1]; $col[10]\n";
        $count++;

        }

    }

}
