#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Extract first N reads (4*N lines) from each *.fastq.gz inside a .tar.gz
# Outputs files named: <member_path_with_slashes_replaced>.head<reads>.fastq

my ($tar_file, $reads, $threads, $help);
$reads   = 10000;  # number of reads (records of 4 lines)
$threads = $ENV{SLURM_CPUS_PER_TASK} // 4;

GetOptions(
    'tar=s'     => \$tar_file,
    'reads=i'   => \$reads,
    'threads=i' => \$threads,
    'help|h'    => \$help,
) or die "Error parsing options\n";

if ($help || !$tar_file) {
    print <<"USAGE";
Usage:
  perl $0 --tar <archive.tar.gz> [--reads 10000] [--threads N]

Description:
  Streams each *.fastq.gz member from the tar archive, decompresses, and writes
  the first <reads> FASTQ records (4*reads lines) to a separate file named:
    <member_path_with_slashes_replaced>.head<reads>.fastq

Options:
  --tar       Path to the .tar.gz archive
  --reads     Number of reads to extract per file (default: 10000)
  --threads   pigz threads (default: SLURM_CPUS_PER_TASK or 4)
  --help      Show this help
USAGE
    exit 1;
}

die "Tar not found: $tar_file\n" unless -e $tar_file;
die "--reads must be >= 1\n" unless $reads >= 1;
die "--threads must be >= 1\n" unless $threads >= 1;

# Convert reads to lines
my $lines = $reads * 4;

# List members once
my @members = `tar -tzf '$tar_file'`;
die "Failed to list archive members for '$tar_file'\n" if $? != 0;
chomp @members;

# Filter for .fastq.gz members
my @fastq_members = grep { /\.fastq\.gz$/ } @members;
die "No .fastq.gz files found in archive '$tar_file'\n" unless @fastq_members;

# Prefer pigz; fall back to zcat
my $have_pigz = system("command -v pigz >/dev/null 2>&1") == 0 ? 1 : 0;
my $decomp_cmd = $have_pigz ? "pigz -p $threads -dc" : "zcat";
warn "Note: pigz not found; using zcat\n" unless $have_pigz;

for my $member (@fastq_members) {
    # Make a safe output filename
    my $safe = $member;
    $safe =~ s#/#_#g;      # replace path separators
    $safe =~ s/\.gz$//;    # drop trailing .gz
    my $out = "${safe}.head${reads}.fastq";

    print "Writing $out from $member ...\n";

    # Stream, decompress, head (4*reads lines), write
    my $cmd = "tar -xOzf '$tar_file' '$member' | $decomp_cmd | head -n $lines > '$out'";
    my $rc = system($cmd);

    if ($rc != 0) {
        warn "Warning: command failed for '$member' (exit $rc). Removing partial '$out'.\n";
        unlink $out;
        next;
    }
}

print "Done. Generated head files for ", scalar(@fastq_members), " member(s).\n";
