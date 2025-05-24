#!/usr/bin/perl
use strict;
use warnings;

# Takes fastq as standard input and return a two column file as standard output
# The output files contains one readID column and one Sequence column
# To use with several .gz files, try zact file1.fq.gz file2.fq.gz ... | extract_reads_fastq.pl
print "readID\tSequence\n"; # Print headers

my $print_next = 0;
my $read_id;

while (<STDIN>) { # Read from standard input
    chomp;
    if (substr($_, 0, 1) eq '@') {
        $read_id = substr($_, 1); # Store the ID without the '@'
        $print_next = 1; # Flag to print the next line (sequence)
    } elsif ($print_next) {
        print "$read_id\t$_\n"; # Print ID and sequence separated by a tab
        $print_next = 0; # Reset the flag
    }
}
