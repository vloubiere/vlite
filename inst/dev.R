setwd("/groups/stark/vloubiere/vlite/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# SLURM wrappers ---------------------------
if(F)
{
  file.edit("R/bsub.R") # bsub gridengine
  file.edit("R/vl_squeue.R")
  file.edit("R/vl_scancel.R")
  file.edit("R/vl_last_err.R")
}

# Dropbox API ------------------------------
if(F)
{
  file.edit("R/dropbox_upload.R")
  file.edit("R/dropbox_download.R")
}

# bamtools ---------------------------------
if(F)
{
  # Import bam
  file.edit("R/importBamRsamtools.R")
  file.edit("R/importBamRaw.R") # Using samtools and fread
}

# Gene ontologies --------------------------
if(F)
{
  file.edit("R/vl_GOenrich.R")
}

# Motifs analyses tools --------------------
if(F)
{
  # Reverse complement DNA sequence
  file.edit("R/revCompDNA.R")

  # pwm converion
  file.edit("R/pwmPercToLog.R")

  # Lasso regression motifs
  file.edit("R/motifLassoRegression.R")

  # Motif counts
  file.edit("R/vl_motifCount.R")
  file.edit("R/vl_motifPos.R")

  # Motif enrichment
  file.edit("R/vl_motifEnrich.R")
  file.edit("R/vl_motifEnrichClusters.R")

  # Download iCisTarget output
  file.edit("R/download_iCistarget.R")
}

# Linear models evaluation ------------------------
if(F)
{
  # Extract model performance stat
  file.edit("R/getModelEquation.R")
  file.edit("R/getModelRMSErsq.R")
  file.edit("R/getModelExpVar.R")
}

# Deeplearning tools -----------------------
if(F)
{
  # Deep learning contrib
  file.edit("R/importContrib.R")
  file.edit("R/contribEnrich.R")
  file.edit("R/contribPlotLogo.R")

  # Modelling diagnostics
  file.edit("R/vl_ROC_AUC.R")
  file.edit("R/vl_PR_AUC.R")
  file.edit("R/vl_PPV.R")
  file.edit("R/vl_mPCC.R")
  file.edit("R/vl_TPR.R")
}

# STRING db tools
if(F)
{
  # Get interactions
  file.edit("R/stringGetDB.R")
  file.edit("R/stringInteraction.R")
  file.edit("R/stringToIgraph.R")

  # Methods to plot
  file.edit("R/plot.vl_STRING.R")
}
# Miscellaneous
if(F)
{
  file.edit("R/gaussianBlur.R")
  file.edit("R/plot_enhancer_motifs.R") # Not reviewed
  file.edit("R/alignSanger.R") # Not reviewed but contains multiple alignment tool
  file.edit("R/select_actMatched_controls.R") # Not reviewed, but contains code usable to sample closest observation
}

# bwtools ---------------------------------
if(F)
{
  # Helper functions
  file.edit("R/helperFunctions_bwTools.R")

  # Coverage
  file.edit("R/bwCoverage.R")
  file.edit("R/bwBinnedCoverage.R")

  # Screenshot
  file.edit("R/bwScreenshot.R")

  # Average tracks
  file.edit("R/bwAverageTrack.R")

  # Bw heatmaps
  file.edit("R/bwHeatmap.R")

}

# BSgenomes tools -------------------------
if(F)
{
  file.edit("R/getBSgenomeSize.R")
  file.edit("R/binBSgenome.R")
  file.edit("R/sampleRegionsBSgenome.R")
  file.edit("R/randomRegionsBSgenome.R")
  file.edit("R/getBSsequence.R")
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
  # Plotting parameters
  file.edit("R/vl_par.R")

  # Function linking values to colors
  file.edit("R/vl_colorRamp.R")

  # Legends
  file.edit("R/tiltAxis.R")
  file.edit("R/vl_legend.R")
  file.edit("R/heatkey.R")

  # Extra labels
  file.edit("R/addRsq.R")
  file.edit("R/addPcc.R")
  file.edit("R/addPval.R")

  # Pie chart
  file.edit("R/vl_pie.R")

  # Boxplot
  file.edit("R/boxplot.R")

  # Heatmap
  file.edit("R/vl_image.R")
  file.edit("R/vl_heatmap.R")

  # Alluvial plot
  file.edit("R/alluvial.R")

  # Balloons plot
  file.edit("R/balloons_key.R")
  file.edit("R/balloons_plot.R")

  # Plot DNA letters/logos
  file.edit("R/plotDNAletter.R") # Used to plot logos
  file.edit("R/addSeqLogo.R")
  file.edit("R/addMotifs.R")

  # MA plot
  file.edit("R/MAplot.R")

  # Methods to plot enrichment analyses output
  file.edit("R/plot_vl_enr.R")
  file.edit("R/plot_vl_enr_clusters.R")

  # Plot table
  file.edit("R/plotTable.R")

  # Scatterplot
  file.edit("R/rasterScatterplot.R")
  file.edit("R/densityScatterplot.R")
  file.edit("R/repelScatterplot.R") # I couldnt completely figure this out yet, careful!

  # Upset plot
  file.edit("R/upsetPlot.R")
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
  file.edit("R/cmd_submit.R")

  # Pipelines (following of wrappers commands)
  file.edit("R/cutnrun_pipeline.R")
  file.edit("R/rnaseq_pipeline.R")
  file.edit("R/orfeome_pipeline.R")
  file.edit("R/orftag_pipeline.R")
  file.edit("R/proseq_pipeline.R")
}
