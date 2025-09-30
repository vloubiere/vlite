#!/usr/bin/perl

# Import required modules
use strict;
use warnings;
use Getopt::Long;

#------------------------------------------------------------------------------
# Global variables for command line options
#------------------------------------------------------------------------------
my $max_reads;
my $help;
my $add_umi;  # Retrieve UMI sequence from i7 and add it to readID_UMI
my $start_seq;  # For PRO-Seq read, only keep reads where R1 starts with the provided eBC sequence
my $trim_length;  # If set with --start_seq, specifies the UMI length that should be extracted after the eBC and appended to the read ID
my $threads;   # Number of threads for pigz

# Parse command line options
GetOptions(
    "n=i"     => \$max_reads,
    "lines=i" => \$max_reads,
    "help|h"  => \$help,
    "umi"     => \$add_umi,
    "start_seq=s" => \$start_seq,
    "trim_length=i" => \$trim_length,
    "threads=i"     => \$threads
) or die "Error in command line arguments\n";

# threads precedence: CLI > SLURM > fallback
$threads = defined($threads)
    ? $threads
    : ($ENV{SLURM_CPUS_PER_TASK} // 4);

die "--threads must be >= 1\n" if $threads < 1;

# Multiply the number of lines by 4 to match fastq formatting
if (defined $max_reads) {
    if ($max_reads < 0) { die "--lines/-n must be non-negative\n"; }
    $max_reads *= 4;  # convert reads to FASTQ lines
}

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
    --start_seq     For PRO-Seq, only keep reads where R1 sequence starts with the provided eBC sequence.
    --trim_length   After trimming --start_seq, extract this many (UMI) nucleotides from R1, append to read ID, and trim from seq/qual.
    -n, --lines     Process only the first N reads from input files.
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
# Build regex patterns for i7, i5 and eBC
#------------------------------------------------------------------------------
my $i7_pattern = join('|', map { quotemeta($_) } @i7_indices); #join indices with a pipe, e.g 'ATCG|TCCA'...
my $i5_pattern = join('|', map { quotemeta($_) } @i5_indices);

# Determine if filtering is needed for i7 and i5
my $use_i7 = $i7_pattern ne '' || $add_umi;  # True if i7_pattern is not empty or --umi is set
my $use_i5 = $i5_pattern ne '';  # True if i5_pattern is not empty

# Determine if filtering is needed for eBC (PRO-Seq)
my $use_start_seq = defined $start_seq && $start_seq ne '';
my $start_seq_re  = $use_start_seq ? qr/^\Q$start_seq\E/ : undef;

# Validate trim_length usage
if (defined $trim_length && $trim_length !~ /^\d+$/) {
    die "Error: --trim_length must be a non-negative integer.\n";
}
if ($trim_length && !$use_start_seq) {
    die "Error: --trim_length requires --start_seq.\n";
}
if ($add_umi && $trim_length) {
    die "Error: --umi and --trim_length cannot be used together.\n";
}

#------------------------------------------------------------------------------
# Verify presence of FASTQ files in tar archive
#------------------------------------------------------------------------------
my @members = `tar -tzf '$tar_file'`;
chomp @members;

# Enforce exactly one R1; others assumed parallel but must exist if required
my @r1_files = grep { /_R1_.*\.fastq\.gz$/ } @members;
die "Error: No R1 FASTQ file found in tar\n" if @r1_files == 0;
die "Error: Multiple R1 FASTQ files found in tar:\n  "
    . join("\n  ", @r1_files)
    . "\nThis script expects exactly ONE R1 file. Please fix the archive or the script."
    if @r1_files > 1;
my ($r1_file) = @r1_files;

# Pick the first matching file for others (must exist if required)
my ($r2_file) = $layout eq 'PAIRED' ? (grep { /_R2_.*\.fastq\.gz$/ } @members) : ();
my ($i1_file) = $use_i7 ? (grep { /_I1_.*\.fastq\.gz$/ } @members) : ();
my ($i2_file) = $use_i5 ? (grep { /_I2_.*\.fastq\.gz$/ } @members) : ();

die "Error: No R2 FASTQ file found in tar (PAIRED layout)\n" if $layout eq 'PAIRED' && !$r2_file;
die "Error: No I1 (i7) FASTQ file found in tar (when --umi or i7 filtering is used)\n" if $use_i7 && !$i1_file;
die "Error: No I2 (i5) FASTQ file found in tar (when i5 filtering is used)\n" if $use_i5 && !$i2_file;

#------------------------------------------------------------------------------
# Open output files
#------------------------------------------------------------------------------
my ($out1, $out2);
if ($layout eq 'PAIRED') {
    open($out1, "|-", "pigz -p $threads -1 > '${fq_prefix}_1.fq.gz'") or die "Cannot open output file 1: $!";
    open($out2, "|-", "pigz -p $threads -1 > '${fq_prefix}_2.fq.gz'") or die "Cannot open output file 2: $!";
} else {
    open($out1, "|-", "pigz -p $threads -1 > '${fq_prefix}.fq.gz'") or die "Cannot open output file: $!";
}

#------------------------------------------------------------------------------
# Open input streams
#------------------------------------------------------------------------------
my $head_cmd = defined($max_reads) ? "| head -n $max_reads" : "";

# Open I1 (i7) stream
my $i1;
if ($use_i7) {
    open($i1, "tar -xOzf '$tar_file' '$i1_file' | pigz -dc $head_cmd |")
      or die "Cannot open I1 file (i7 reads): $!\n";
}

# Open I2 (i5) stream if needed
my $i2;
if ($use_i5) {
    open($i2, "tar -xOzf '$tar_file' '$i2_file' | pigz -dc $head_cmd |")
        or die "Cannot open I2 file (i5 reads): $!\n";
}

# Open R1 stream
open(my $r1, "tar -xOzf '$tar_file' '$r1_file' | pigz -dc $head_cmd |")
    or die "Cannot open R1 file: $!\n";

# Open R2 stream if paired-end
my $r2;
if ($layout eq 'PAIRED') {
    open($r2, "tar -xOzf '$tar_file' '$r2_file' | pigz -dc $head_cmd |")
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

    $total_reads++; # One read has been successfully imported

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
    my $r1_seq_ok = 1; # Default to pass if no eBC start_seq provided
    if ($use_start_seq) {
        my $r1_seq = $r1_chunk[1]; chomp $r1_seq;
        $r1_seq_ok = ($r1_seq =~ $start_seq_re) ? 1 : 0;
    }

    # Compute PRO-Seq eBC length if provided
    my $eBC_len = 0;
    if ($use_start_seq) {
        my $r1_seq = $r1_chunk[1]; chomp $r1_seq;
        if ($r1_seq =~ $start_seq_re) {
            $eBC_len = length($start_seq);
        }
    }

    # If both indices match (and eBC if start_seq provided) match, write the reads
    if ($i7_ok && $i5_ok && $r1_seq_ok) {
        $matched_reads++; # Increment matching read counts

        # Append i7 to read ID if --umi is set
        if ($add_umi) {
            $r1_chunk[0] =~ s/^(\S+)/$1_$current_i7_seq/;  # Append i7 (UMI) to the first block of R1 ID
            if ($layout eq 'PAIRED') {
                $r2_chunk[0] =~ s/^(\S+)/$1_$current_i7_seq/;  # Append i7 (UMI) to the first block of R2 ID
            }
        }

        # If start_seq provided (PRO-Seq eBC), trim it from R1 sequence and quality
        if ($eBC_len > 0) {
            # Trim sequence (R1 line 2)
            my $seq = $r1_chunk[1]; chomp $seq;
            substr($seq, 0, $eBC_len, '');
            $r1_chunk[1] = "$seq\n";

            # Trim quality (R1 line 4)
            my $qual = $r1_chunk[3]; chomp $qual;
            substr($qual, 0, $eBC_len, '');
            $r1_chunk[3] = "$qual\n";
        }

        # If --trim_length is set, extract UMI from the start of (post-eBC) R1 and append to IDs
        if ($use_start_seq && $trim_length && $trim_length > 0) {
            # Extract UMI from R1 sequence and trim sequence/quality accordingly
            my $seq  = $r1_chunk[1]; chomp $seq;
            my $qual = $r1_chunk[3]; chomp $qual;

            # Guard against short reads
            my $take = $trim_length;
            $take = length($seq) if $take > length($seq);
            my $umi = substr($seq, 0, $take);

            substr($seq,  0, $take, '');
            substr($qual, 0, $take, '');

            $r1_chunk[1] = "$seq\n";
            $r1_chunk[3] = "$qual\n";

            # Append UMI to first token of read IDs
            $r1_chunk[0] =~ s/^(\S+)/$1_$umi/;
            if ($layout eq 'PAIRED') {
                $r2_chunk[0] =~ s/^(\S+)/$1_$umi/;
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
if ($layout eq 'PAIRED') {
    close($r2);
}

close($out1) or die "pigz output 1 failed: $!";
if ($layout eq 'PAIRED') {
    close($out2) or die "pigz output 2 failed: $!";
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
print "R1 start filter: ", ($use_start_seq ? $start_seq : "Disabled"), "\n";
print "Trim-length UMI: ", ($use_start_seq && $trim_length ? $trim_length : "Disabled"), "\n";
print "Total reads: $total_reads\n";
print "Matched reads: $matched_reads\n";
printf "Match rate: %.2f%%\n",
    ($total_reads > 0 ? ($matched_reads / $total_reads * 100) : 0);
