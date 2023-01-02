#!/usr/bin/perl

## Using the strict + warning pragmas to enforce good programming practices
use warnings;   ## Useful to find potential problems in the code
use strict;     ## Requires the use of properly defined lexical and package variables

## Writing a short readme:
## $0 is a special variable in Perl that stores the name of the script being run
my $usage = "$0 input_file.fastq output_file.fasta\n";  ## $0 is a special variable in perl
die "\n$usage\n" unless @ARGV;

## Working with filehandles
open FASTQ, "<", $ARGV[0]; ## Input file in FASTQ format
open FASTA, "<", $ARGV[1]; ## Output file (we will create it in FASTA format)

## Notes on the FASTQ input file:
## Sequences in FASTQ format are stored across 4 lines:
# line 1 = sequence name (starts with an @)
# line 2 = sequence
# line 3 = spacer (starts with a +; not used in FASTA format)
# line 4 = quality scores (not used in FASTA format)

## Creating variables to be used later
my $seq_name;
my $sequence;
my $line_counter;

## Iterating through the file line per line
while (my $line = <IN>){

	chomp $lines;
	$line_counter++;

	## Using modulo to keep track of which line we are looking at:
	my $modulo = $line_counter % 4;

	## modulo = 1 => line with sequence name
	## modulo = 2 => line with sequence
	## modulo = 3 => line with spacer (+)
	## modulo = 0 => line with quality score; 4 % 4 = 0, no remainder to the division

	if ($modulo == 1){
		$seq_name = $line;
	}

	elsif ($modulo == 2){

		$sequence = $line;

		## Breaking string to 60 ASCII characters per line with unpack()
		## Sequence lines in FASTA format are usually 60 characters-wide
		my @seq60_wide = unpack ("(A60)*", $line);

		## Printing FASTA header
		print OUT ">$seq_name\n";

		## Printing FASTA sequence
		while (my $seq60 = shift@seq50_wide){
			print FASTA "$seq60\n";
		}

	}

}