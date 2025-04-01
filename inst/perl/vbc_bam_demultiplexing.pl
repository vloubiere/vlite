#!/usr/bin/perl

# Script: vbc_bam_demultiplexing.pl
# Description: Process and filter sequencing reads based on index sequences and read properties
# Version: 1.0.3 (January 2024)

use strict;
use warnings;
use IO::Compress::Gzip qw(gzip $GzipError);

# Print help if --help or -h is provided
if (@ARGV == 0 || grep {$_ eq '--help' || $_ eq '-h'} @ARGV) {
    print_help();
    exit;
}

# Subroutine to print comprehensive help documentation
sub print_help {
    print <<'END_HELP';
NAME
    vbc_bam_demultiplexing.pl - Process and filter sequencing reads based on index sequences and read properties

SYNOPSIS
    ./vbc_bam_demultiplexing.pl <PAIRED|SINGLE> <i7 index(es)> <i5 index(es)> <i7_column> <i5_column> <output_prefix> [start_seq] [trim_length]

DESCRIPTION
    This script processes SAM/BAM format sequencing data from standard input, filtering reads based on
    index sequences (i7 and i5) and optionally processing read sequences for specific protocols like PROseq.
    It can handle both paired-end and single-end data.
    For paired-end data, index information is expected only in the first read of each pair.

REQUIRED ARGUMENTS
    <PAIRED|SINGLE>    Specify the read type:
                       PAIRED - for paired-end sequencing data
                       SINGLE - for single-end sequencing data

    <i7 index(es)>     The i7 index sequence(s) to match. Multiple indexes should be comma-separated (this is only useful to analyze ORFtag).
                       Use 'none' if no i7 index filtering is needed.
                       Example: "AAGGCCTT" or "AAGGCCTT,CCTTAAGG" or "none"
                       Note: BC:Z: prefix is automatically added

    <i5 index(es)>     The i5 index sequence(s) to match. Multiple indexes should be comma-separated (this is only useful to analyze ORFtag).
                       Use 'none' if no i5 index filtering is needed.
                       Example: "ACGTACGT" or "ACGTACGT,TGCATGCA" or "none"
                       Note: B2:Z: prefix is automatically added

    <i7_column>        1-based column number in the SAM file containing the i7 index
                       Typically 14 in standard format

    <i5_column>        1-based column number in the SAM file containing the i5 index
                       Typically 12 in standard format

    <output_prefix>    Prefix for output files. The script will append:
                       '_1.fq.gz' and '_2.fq.gz' for paired-end data
                       '.fq.gz' for single-end data

OPTIONAL ARGUMENTS
    [start_seq]        Sequence that must be present at the start of reads
                       If provided, only reads starting with this sequence will be kept
                       The sequence will be trimmed from matching reads

    [trim_length]      Number of nucleotides to extract after trimming start_seq
                       Only applied if start_seq is provided
                       The extracted sequence will be added to the read ID
                       Remaining sequence will be kept as the read

END_HELP
}

# Validate command line arguments
if (@ARGV < 6) {
    die "vbc_bam_demultiplexing.pl: Insufficient arguments.\nUse --help or -h for usage information.\n";
}

# Get and validate read type (PAIRED or SINGLE)
my $read_type = uc($ARGV[0]);
unless ($read_type eq 'PAIRED' || $read_type eq 'SINGLE') {
    die "First argument must be either PAIRED or SINGLE\n";
}

# Process i7 indexes (now first)
my @i7_indexes = ();
if ($ARGV[1] && $ARGV[1] ne 'none') {
    @i7_indexes = split(',', $ARGV[1]);
    @i7_indexes = map { "^BC:Z:" . $_ } @i7_indexes;
}

# Process i5 indexes (now second)
my $i5_input = $ARGV[2];
my @i5_indexes = ();
if ($i5_input && $i5_input ne 'none') {
    @i5_indexes = split(',', $i5_input);
    @i5_indexes = map { "^B2:Z:" . $_ } @i5_indexes;
}

# Convert column numbers to 0-based indexing for internal use
my $i7_column = $ARGV[3] - 1;
my $i5_column = $ARGV[4] - 1;
my $output_prefix = $ARGV[5];

# Process optional arguments
my $start_seq = $ARGV[6] if @ARGV > 6;
my $trim_length = $ARGV[7] if @ARGV > 7;

# Validate inputs
if ($i7_column < 0 || $i5_column < 0) {
    die "Column numbers must be positive integers\n";
}

# Open output file(s) based on read type
my ($out1, $out2);
if ($read_type eq 'PAIRED') {
    $out1 = IO::Compress::Gzip->new("${output_prefix}_1.fq.gz")
        or die "Cannot open output file 1: $GzipError";
    $out2 = IO::Compress::Gzip->new("${output_prefix}_2.fq.gz")
        or die "Cannot open output file 2: $GzipError";
} else {
    $out1 = IO::Compress::Gzip->new("${output_prefix}.fq.gz")
        or die "Cannot open output file: $GzipError";
}

# Track if any reads are found
my $reads_found = 0;

# Process input SAM/BAM file from STDIN
while (my $line = <STDIN>) {
    chomp $line;

    # Skip header lines
    next if ($line =~ /^\@/);

    my @fields = split("\t", $line);

    # Determine if we should process this read based on read type
    my $process_read = 0;
    if ($read_type eq 'PAIRED') {
        $process_read = ($fields[1] & 0x40);  # Check if first mate in pair
    } else {
        $process_read = !($fields[1] & 0x1);  # Check if not paired
    }

    if ($process_read) {
        # Initialize match flag
        my $matches = 1;

        # Only check indexes on first mate in paired-end mode or all reads in single-end mode
        my $should_check_indexes = ($read_type eq 'PAIRED') ? ($fields[1] & 0x40) : 1;

        if ($should_check_indexes) {
            # Validate column existence for index checking
            if (@i7_indexes > 0 && @fields <= $i7_column) {
                die "Error: i7 column (${\($i7_column + 1)}) exceeds number of fields in input. Line: $.\n";
            }
            if (@i5_indexes > 0 && @fields <= $i5_column) {
                die "Error: i5 column (${\($i5_column + 1)}) exceeds number of fields in input. Line: $.\n";
            }

            # Check i7 indexes first if provided
            if (@i7_indexes > 0) {
                my $i7_match = 0;
                foreach my $i7_index (@i7_indexes) {
                    if ($fields[$i7_column] =~ /$i7_index/) {
                        $i7_match = 1;
                        last;
                    }
                }
                $matches = $i7_match;
            }

            # Check i5 indexes if i7 matched and i5 indexes are provided
            if ($matches && @i5_indexes > 0) {
                my $i5_match = 0;
                foreach my $i5_index (@i5_indexes) {
                    if ($fields[$i5_column] =~ /$i5_index/) {
                        $i5_match = 1;
                        last;
                    }
                }
                $matches = $i5_match;
            }
        }

        # Check start sequence if provided
        if ($matches && $start_seq) {
            $matches = ($fields[9] =~ /^$start_seq/);
        }

        if ($matches) {
            $reads_found = 1;

            # Process read sequence and quality scores
            my $seq = $fields[9];
            my $qual = $fields[10];
            my $read_id = $fields[0];

            # Handle sequence processing if start_seq and trim_length are provided
            if ($start_seq && $trim_length) {
                # First trim the start sequence if it matches
                if ($seq =~ /^$start_seq/) {
                    # Remove the start sequence
                    $seq =~ s/^$start_seq//;
                    $qual = substr($qual, length($start_seq));

                    # Extract the next trim_length nucleotides
                    my $extracted_seq = substr($seq, 0, $trim_length);

                    # Keep the remaining sequence (after trim_length nucleotides)
                    $seq = substr($seq, $trim_length);
                    $qual = substr($qual, $trim_length);

                    # Append the extracted sequence to the read ID
                    $read_id .= ":" . $extracted_seq;
                }
            }

            # Output reads based on read type
            if ($read_type eq 'PAIRED') {
                # Write first mate
                print $out1 "\@$read_id\n$seq\n+\n$qual\n";

                # Process and write second mate
                my $second_mate = <STDIN>;
                chomp $second_mate;
                my @second_fields = split("\t", $second_mate);
                print $out2 "\@" . $second_fields[0] . "\n";
                print $out2 $second_fields[9] . "\n";
                print $out2 "+\n";
                print $out2 $second_fields[10] . "\n";
            } else {
                # Write single-end read
                print $out1 "\@$read_id\n$seq\n+\n$qual\n";
            }
        }
    }
}

# Close output files
close($out1);
close($out2) if $read_type eq 'PAIRED';

# Check if any reads were found and provide appropriate error message if not
if (!$reads_found) {
    my $error_message = $read_type eq 'PAIRED' ?
        "No paired-end reads detected in the SAM file." :
        "No single-end reads detected in the SAM file.";
    die $error_message;
}

# Print completion message with appropriate file names
print "Processing complete. Output file(s): ";
print $read_type eq 'PAIRED' ?
    "${output_prefix}_1.fq.gz and ${output_prefix}_2.fq.gz\n" :
    "${output_prefix}.fq.gz\n";
