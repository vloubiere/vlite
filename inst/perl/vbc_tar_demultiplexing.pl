#!/usr/bin/perl

# Import required modules
use strict;
use warnings;
use IO::Compress::Gzip qw(gzip $GzipError);
use Getopt::Long;

#------------------------------------------------------------------------------
# Global variables for command line options
#------------------------------------------------------------------------------
my $max_reads;
my $help;
my $add_umi;  # Retrieve UMI sequence from i7 and add it to readID_UMI

# Parse command line options
GetOptions(
    "n=i"     => \$max_reads,
    "lines=i" => \$max_reads,
    "help|h"  => \$help,
    "umi"     => \$add_umi
) or die "Error in command line arguments\n";

# Display help message if requested or insufficient arguments
if ($help || @ARGV < 4) {
    print <<'END_HELP';
Usage: ./script.pl [options] <tar_file> <i7_indices> <i5_indices> <output_prefix> [layout]

Description:
    This script processes FASTQ files from a tar archive, filtering reads based on
    specified i7 and i5 index sequences. It can handle both single-end and paired-end reads, as well as i7 UMIs.

    Index sequences are read from separate FASTQ files:
    - *_I1_*.fastq.gz for i7 sequences
    - *_I2_*.fastq.gz for i5 sequences
    - *_R1_*.fastq.gz for read 1
    - *_R2_*.fastq.gz for read 2 (paired-end only)

Arguments:
    tar_file        Path to the tar.gz archive containing FASTQ files
    i7_indices      Comma-separated list of i7 index sequences (use 'none' if not applicable).
    i5_indices      Comma-separated list of i5 index sequences (use 'none' if not applicable).
    output_prefix   Prefix for output files
    layout          [Optional] Read layout: 'SINGLE' or 'PAIRED' (default: PAIRED)

Options:
    --umi           The UMI is extracted from the i7 read (I1) and appended to the read ID,
                    separated by an underscore (_).
    -n, --lines     Process only the first N lines from input files
                    Note: Each read consists of 4 lines, so N=1000 processes 250 reads
    -h, --help      Display this help message and exit

Output:
    For paired-end data:
        <output_prefix>_1.fq.gz
        <output_prefix>_2.fq.gz
    For single-end data:
        <output_prefix>.fq.gz

Example:
    ./script.pl -n 1000 input.tar.gz ATTCAGAA,CTGAAGCT none output_files PAIRED

Notes:
    - The script expects four gzipped FASTQ files in the tar archive
    - Index sequences are case-sensitive
END_HELP
    exit;
}

#------------------------------------------------------------------------------
# Parse command line arguments
#------------------------------------------------------------------------------
my ($tar_file, $i7_list, $i5_list, $fq_prefix, $layout) = @ARGV;
die "Usage: $0 <tar_file> <i7_list> <i5_list> <fq_prefix> <layout>\n" unless defined $layout;

# Input validation
die "Error: Tar file '$tar_file' not found\n" if !-e $tar_file;
if (lc($i7_list) eq 'none' && lc($i5_list) eq 'none') {
    die "Error: i7 and i5 are both set to 'none'.\n";
}
$layout = uc($layout); # Layout to upper case
die "Error: Invalid layout '$layout'\n" if $layout ne 'PAIRED' && $layout ne 'SINGLE';

# Check if i7 or i5 is 'none' and set them to empty strings
if (lc($i7_list) eq 'none') {
    $i7_list = '';
}
if (lc($i5_list) eq 'none') {
    $i5_list = '';
}

# Validate --umi usage
if ($add_umi) {
    if ($i7_list ne '') {
        die "Error: --umi option requires i7 to be set to 'none' or an empty string.\n";
    }
}

# If several comma-separated indexes were provided, split them
my @i7_indices = split(/,/, $i7_list);
my @i5_indices = split(/,/, $i5_list);

#------------------------------------------------------------------------------
# Build regex patterns for i7 and i5
#------------------------------------------------------------------------------
my $i7_pattern = join('|', map { quotemeta($_) } @i7_indices); #join indices with a pipe, e.g 'ATCG|TCCA'...
my $i5_pattern = join('|', map { quotemeta($_) } @i5_indices);

# Determine if filtering is needed for i7 and i5
my $use_i7 = $i7_pattern ne '' || $add_umi;  # True if i7_pattern is not empty or --umi is set
my $use_i5 = $i5_pattern ne '';  # True if i5_pattern is not empty

#------------------------------------------------------------------------------
# Verify presence of FASTQ files in tar archive
#------------------------------------------------------------------------------
my $check_r1 = `tar -tzf $tar_file | grep '_R1_.*fastq.gz\$'`;
die "Error: No R1 FASTQ files found\n" if !$check_r1;

if ($layout eq 'PAIRED') {
    my $check_r2 = `tar -tzf $tar_file | grep '_R2_.*fastq.gz\$'`;
    die "Error: No R2 FASTQ files found\n" if !$check_r2;
}

if ($use_i7) {
    my $check_i1 = `tar -tzf $tar_file | grep '_I1_.*fastq.gz\$'`;
    die "Error: No I1 FASTQ files (i7 sequences) found\n" if !$check_i1;
}

if ($use_i5) {
    my $check_i2 = `tar -tzf $tar_file | grep '_I2_.*fastq.gz\$'`;
    die "Error: No I2 FASTQ files (i5 sequences) found\n" if !$check_i2;
}

#------------------------------------------------------------------------------
# Open output files
#------------------------------------------------------------------------------
my ($out1, $out2);
if ($layout eq 'PAIRED') {
    $out1 = IO::Compress::Gzip->new("${fq_prefix}_1.fq.gz")
        or die "Cannot open output file 1: $GzipError\n";
    $out2 = IO::Compress::Gzip->new("${fq_prefix}_2.fq.gz")
        or die "Cannot open output file 2: $GzipError\n";
} else {
    $out1 = IO::Compress::Gzip->new("${fq_prefix}.fq.gz")
        or die "Cannot open output file: $GzipError\n";
}

#------------------------------------------------------------------------------
# Open input streams
#------------------------------------------------------------------------------
my $head_cmd = defined($max_reads) ? "| head -n $max_reads" : "";

# Open I1 (i7) stream
my $i1;
if ($use_i7) {
    open($i1, "tar -xOzf $tar_file --wildcards '*_I1_*fastq.gz' | zcat $head_cmd |")
      or die "Cannot open I1 file (i7 reads): $!\n";
}

# Open I2 (i5) stream if needed
my $i2;
if ($use_i5) {
    open($i2, "tar -xOzf $tar_file --wildcards '*_I2_*fastq.gz' | zcat $head_cmd |")
        or die "Cannot open I2 file (i5 reads): $!\n";
}

# Open R1 stream
open(my $r1, "tar -xOzf $tar_file --wildcards '*_R1_*fastq.gz' | zcat $head_cmd |")
    or die "Cannot open R1 file: $!\n";

# Open R2 stream if paired-end
my $r2;
if ($layout eq 'PAIRED') {
    open($r2, "tar -xOzf $tar_file --wildcards '*_R2_*fastq.gz' | zcat $head_cmd |")
        or die "Cannot open R2 file: $!\n";
}

#------------------------------------------------------------------------------
# Main processing loop
#------------------------------------------------------------------------------
my ($total_reads, $matched_reads) = (0, 0);

while (1) {
    # Read I1 chunk (i7)
    my @i1_chunk = ();
    if ($use_i7) {
        for (1..4) {
            my $line = <$i1>;
            last unless defined $line;
            push @i1_chunk, $line;
        }
        last if scalar(@i1_chunk) < 4;
    }

    # Read I2 chunk (i5) if needed
    my @i2_chunk = ();
    if ($use_i5) {
        for (1..4) {
            my $line = <$i2>;
            last unless defined $line;
            push @i2_chunk, $line;
        }
        last if scalar(@i2_chunk) < 4;
    }

    # Read R1 chunk
    my @r1_chunk = ();
    for (1..4) {
        my $line = <$r1>;
        last unless defined $line;
        push @r1_chunk, $line;
    }
    last if scalar(@r1_chunk) < 4;

    # Read R2 chunk if paired-end
    my @r2_chunk = ();
    if ($layout eq 'PAIRED') {
        for (1..4) {
            my $line = <$r2>;
            last unless defined $line;
            push @r2_chunk, $line;
        }
        last if scalar(@r2_chunk) < 4;
    }

    $total_reads++; # One read has been succesfully imported

    # Extract sequences from I1 and I2
    my $current_i7_seq = '';
    if ($use_i7) {
        chomp($current_i7_seq = $i1_chunk[1]); # chomp removes any trailing newline
    }
    my $current_i5_seq = '';
    if ($use_i5) {
        chomp($current_i5_seq = $i2_chunk[1]);
    }

    # Check for matches at the start of reads
    my $i7_ok = 1;  # Default to pass if not using i7
    if ($use_i7 && !$add_umi) { # If i7 is a UMI, this regexpr matching is skipped here
        $i7_ok = ($current_i7_seq =~ /^(?:$i7_pattern)/) ? 1 : 0;
    }
    my $i5_ok = 1;  # Default to pass if not using i5
    if ($use_i5) {
        $i5_ok = ($current_i5_seq =~ /^(?:$i5_pattern)/) ? 1 : 0;
    }

    # If both indices match, write the reads
    if ($i7_ok && $i5_ok) {
        $matched_reads++; # Increment matching read counts

        # Append i7 to read ID if --umi is set
        if ($add_umi) {
            $r1_chunk[0] =~ s/^(\S+)/$1_$current_i7_seq/;  # Append i7 (UMI) to the first block of R1 ID
            if ($layout eq 'PAIRED') {
                $r2_chunk[0] =~ s/^(\S+)/$1_$current_i7_seq/;  # Append i7 (UMI) to the first block of R2 ID
            }
        }

        # Write reads to output files
        if ($layout eq 'SINGLE') {
            print $out1 @r1_chunk;
        } else {
            print $out1 @r1_chunk;
            print $out2 @r2_chunk;
        }
    }
}

#------------------------------------------------------------------------------
# Close file handles
#------------------------------------------------------------------------------
if ($use_i7) { close($i1); }
if ($use_i5) { close($i2); }

close($r1);
close($out1);
if ($layout eq 'PAIRED') {
    close($r2);
    close($out2);
}

#------------------------------------------------------------------------------
# Print summary
#------------------------------------------------------------------------------
print "Processing complete.\n";
print "Files: ",
    ($layout eq 'PAIRED'
        ? "${fq_prefix}_1.fq.gz, ${fq_prefix}_2.fq.gz"
        : "${fq_prefix}.fq.gz"),
    "\n";
print "UMI mode: ", ($add_umi ? "Enabled" : "Disabled"), "\n";
print "Total reads: $total_reads\n";
print "Matched reads: $matched_reads\n";
printf "Match rate: %.2f%%\n",
    ($total_reads > 0 ? ($matched_reads / $total_reads * 100) : 0);
