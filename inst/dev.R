setwd("/groups/stark/vloubiere/vlite/")

# bwtools ---------------------------------
if(F)
{
  # Helper functions
  file.edit("R/helperFunctions_bwTools.R")

  # Coverage
  file.edit("R/bwCoverage.R")
  file.edit("R/bwBinnedCoverage.R")
  file.edit("R/bwAverageTrack.R")

  # Screenshot
  file.edit("R/bwScreenshot.R")
}

# bedtools --------------------------------
if(F)
{
  # Helper functions
  file.edit("R/helperFunctions_bedTools.R")

  # Import/export
  file.edit("R/importBed.R")
  file.edit("R/exportBed.R")

  # Standalone funcitons
  file.edit("R/resizeBed.R")
  file.edit("R/collapseBed.R")
  file.edit("R/closestBed.R")
  file.edit("R/binBed.R")

  # Fuinctions using overlapBed
  file.edit("R/overlapBed.R")
  file.edit("R/covBed.R")
  file.edit("R/intersectBed.R")
  file.edit("R/subtractBed.R")
  file.edit("R/clipBed.R")
}

# plots --------------------------------
if(F)
{
  file.edit("R/gPar.R")
  file.edit("R/tiltAxis.R")
  file.edit("R/gLegend.R")
  file.edit("R/gHeatkey.R")

  file.edit("R/gImage.R")
  file.edit("R/gHeatmap.R")
}

# Genomics pipelines ---------------------
if(F)
{
  # Helper functions
  file.edit("R/helperFunctions_genomicPipelines.R")

  # Donwload data
  # SRA toolkit
  file.edit("R/cmd_downloadSRA.R")
  # fqs (or other...)
  file.edit("R/cmd_download.R")

  # Demultiplexing
  file.edit("R/cmd_demultiplexVBCfile.R") # Wrapper
  file.edit("inst/perl/vbc_tar_demultiplexing.pl") # Perl subscript
  file.edit("inst/perl/vbc_bam_demultiplexing.pl") # Perl subscript

  # Trimming
  # Illumina adaptors
  file.edit("R/cmd_trimIlluminaAdaptors.R")
  # Custom adaptors (PROseq)
  file.edit("R/cmd_trimProseqAdaptors.R")

  # Alignment
  # Bowtie1 alignment
  file.edit("R/cmd_alignBowtie.R")
  # Bowtie2 alignment
  file.edit("R/cmd_alignBowtie2.R")
  # Rsubread alignment (RNA-Seq)
  file.edit("R/cmd_alignRsubread.R") # Wrapper
  file.edit("inst/Rscripts/align_Rsubread.R") # R subscript

  # Post-alignment processing
  # Collapse bam file (ORFtag)
  file.edit("R/cmd_collapseBam.R")
  # Extract unaligned reads (PRO-Seq spike-in)
  file.edit("R/cmd_extractUnalignedReads.R")

  # Count reads
  # Rsubread (RNA-Seq)
  file.edit("R/cmd_countRsubread.R") # Wrapper
  file.edit("inst/Rscripts/count_Rsubread.R") # R subscript
  # Per Barcode (ORFeome)
  file.edit("R/cmd_countBCreads.R") # Wrapper
  file.edit("inst/Rscripts/BC_counts.R") # R subscript
  # Assign reads to exons (ORFtag)
  file.edit("R/cmd_assignInsertions.R") # Wrapper
  file.edit("inst/Rscripts/assign_ORFtag_insertions.R") # R subscript
  # UMI collapsing (from bam)
  file.edit("R/cmd_umiCountsFromBam.R") # R subscript
  file.edit("inst/Rscripts/umiCountsFromBam.R") # R subscript
  # UMI counts (PROseq)
  file.edit("R/cmd_umiCountsProseq.R")
  file.edit("inst/Rscripts/umiCountsProseq.R")

  # Peak calling (ChIP-Seq)
  file.edit("R/cmd_peakCalling.R") # Wrapper
  file.edit("R/cmd_confidentPeaks.R") # Wrapper
  file.edit("inst/Rscripts/confident_peaks.R") # R subscript

  # Generate bigwig tracks
  # Bedgraph to bigwig
  file.edit("R/cmd_bedgraphToBigwig.R") # Wrapper
  file.edit("inst/Rscripts/bedgraph_to_bigwig.R") # R subscript
  # BAM to bigwig
  file.edit("R/cmd_bamToBigwig.R") # Wrapper
  file.edit("inst/Rscripts/bam_to_bigwig.R") # R subscript
  # BED to bigwig
  file.edit("R/cmd_bedToBigwig.R") # Wrapper
  file.edit("inst/Rscripts/bed_to_bigwig.R") # R subscript
  # PROSeq UMI counts to bigwig
  file.edit("R/cmd_umiToBigwigProseq.R") # Wrapper

  # Submit commands
  file.edit("R/cmd_submit.R") # Wrapper
  file.edit("R/bsub.R") # bsub gridengine

  # Pipelines (following of wrappers commands)
  file.edit("R/cutnrun_pipeline.R")
  file.edit("R/rnaseq_pipeline.R")
  file.edit("R/orfeome_pipeline.R")
  file.edit("R/orftag_pipeline.R")
  file.edit("R/proseq_pipeline.R")
}
