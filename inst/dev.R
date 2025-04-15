setwd("/groups/stark/vloubiere/vlite/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
# file.edit("inst/test.R")

# Tests/wip --------------------------------
if(F)
{
  file.edit("inst/test/newUmiCollapsingStrat.R")
}

# bamtools ---------------------------------
if(F)
{
  # Import bam
  file.edit("R/importBamRsamtools.R")
  file.edit("R/importBamRaw.R") # Using samtools and fread
}

# bedtools ---------------------------------
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

  # Functions using overlapBed
  file.edit("R/overlapBed.R")
  file.edit("R/covBed.R")
  file.edit("R/intersectBed.R")
  file.edit("R/subtractBed.R")
  file.edit("R/clipBed.R")
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

# Motifs analyses tools --------------------
if(F)
{
  # Reverse complement DNA sequence
  file.edit("R/revCompDNA.R")

  # pwm converion
  file.edit("R/pwmPercToLog.R")

  # Motif counts
  file.edit("R/vl_motifCount.R")
  file.edit("R/vl_motifPos.R")

  # Motif enrichment
  file.edit("R/vl_motifEnrich.R")

  # Lasso regression motifs
  file.edit("R/motifLassoRegression.R")

  # Download iCisTarget output
  file.edit("R/download_iCistarget.R")
}

# Gene ontologies --------------------------
if(F)
{
  file.edit("R/vl_GOenrich.R")
}

# SLURM wrappers ---------------------------
if(F)
{
  file.edit("R/bsub.R") # bsub gridengine
  file.edit("R/vl_submit.R")
  file.edit("R/vl_squeue.R")
  file.edit("R/vl_scancel.R")
  file.edit("R/vl_last_err.R")
}

# Download data ----------------------------
if(F)
{
  # SRA toolkit
  file.edit("R/cmd_downloadSRA.R")
  # fqs (or other...)
  file.edit("R/cmd_download.R")
}

# Dropbox API ------------------------------
if(F)
{
  file.edit("R/dropboxUpload.R")
  file.edit("R/dropboxDownload.R")
}

# Demultiplexing commands ------------------
if(F)
{
  file.edit("R/cmd_demultiplexVBCfile.R") # Wrapper
  file.edit("inst/perl/vbc_tar_demultiplexing.pl") # Perl subscript
  file.edit("inst/perl/vbc_bam_demultiplexing.pl") # Perl subscript
}

# Pipelines commands -----------------------
if(F)
{
  # Trimming ----
  # Illumina adaptors
  file.edit("R/cmd_trimIlluminaAdaptors.R")
  # Custom adaptors (PROseq)
  file.edit("R/cmd_trimProseqAdaptors.R")

  # Alignment ----
  # Bowtie1 alignment
  file.edit("R/cmd_alignBowtie.R")
  # Bowtie2 alignment
  file.edit("R/cmd_alignBowtie2.R")
  # Rsubread alignment (RNA-Seq)
  file.edit("R/cmd_alignRnaRsubread.R") # Wrapper
  file.edit("inst/Rscript/align_rna_Rsubread.R") # R subscript

  # BAM post-processing ----
  # Collapse bam file unique mapping positions (ORFtag)
  file.edit("R/cmd_collapseBam.R")
  # UMI collapsing from bam (human STARR-Seq, ...)
  file.edit("R/cmd_umiCountsFromBam.R") # Wrapper
  file.edit("inst/Rscript/umiCountsFromBam.R") # R subscript
  # Extract unaligned reads (PRO-Seq spike-in)
  file.edit("R/cmd_extractUnalignedReadsFromBam.R")
  # UMI counts from PROseq bam
  file.edit("R/cmd_umiCollapsingProseq.R") # Wrapper
  file.edit("inst/Rscript/umiCollapsingProseq.R") #R subscript

  # Count reads per feature ----
  # Rsubread (RNA-Seq)
  file.edit("R/cmd_countRsubread.R") # Wrapper
  file.edit("inst/Rscript/count_Rsubread.R") # R subscript
  # Count ORFeomes barcode (ORFeome)
  file.edit("R/cmd_countBCreads.R") # Wrapper
  file.edit("inst/Rscript/BC_counts.R") # R subscript
  # Assign reads to exons (ORFtag)
  file.edit("R/cmd_assignInsertions.R") # Wrapper
  file.edit("inst/Rscript/assign_ORFtag_insertions.R") # R subscript
  # From .rds annotations (PRO-Seq)
  file.edit("R/cmd_countPROseqReads.R") # Wrapper
  file.edit("inst/Rscript/count_PROseq_reads.R") # R subscript

  # Generate bigwigs ----
  # Bedgraph to bigwig
  file.edit("R/cmd_bedgraphToBigwig.R") # Wrapper
  file.edit("inst/Rscript/bedgraph_to_bigwig.R") # R subscript
  # BAM to bigwig
  file.edit("R/cmd_bamToBigwig.R") # Wrapper
  file.edit("inst/Rscript/bam_to_bigwig.R") # R subscript
  # BED to bigwig
  file.edit("R/cmd_bedToBigwig.R") # Wrapper
  file.edit("inst/Rscript/bed_to_bigwig.R") # R subscript
  # PROSeq UMI counts to bigwig
  file.edit("R/cmd_umiToBigwigProseq.R") # Wrapper
  file.edit("inst/Rscript/umiToBigwigProseq.R") # R subscript
}

# Genomics pipelines -----------------------
if(F)
{
  # Helper functions
  file.edit("R/helperFunctions_genomicPipelines.R")

  # Pipelines ----
  file.edit("R/cutnrunProcessing.R")
  file.edit("R/rnaseqProcessing.R")
  file.edit("R/orfeomeProcessing.R")
  file.edit("R/orftagProcessing.R")
  file.edit("R/proseqProcessing.R")
}

# Peak calling commands --------------------
if(F)
{
  file.edit("R/cmd_peakCalling.R") # Wrapper
  file.edit("R/cmd_confidentPeaks.R") # Wrapper
  file.edit("inst/Rscript/confident_peaks.R") # R subscript
}

# DESeq2 commands --------------------------
if(F)
{
  # RNA-Seq
  file.edit("R/cmd_DESeq2.R") # Wrapper
  file.edit("inst/Rscript/DESeq2_analysis.R") # R subscript
  # PRO-Seq (spike-in norm...)
  file.edit("R/cmd_DESeq2_PROseq.R") # Wrapper
  file.edit("inst/Rscript/DESeq2_PROseq_analysis.R") # R subscruot
}

# MAGECK commands --------------------------
if(F)
{
  file.edit("R/cmd_MAGECK_ORFeome.R") # Wrapper
  file.edit("inst/Rscript/compute_MAGECK_count_tables_ORFeome.R") # R subscript
  file.edit("inst/Rscript/volcano_plots_MAgECK.R") # R subscript
  file.edit("inst/Rscript/merge_gene_summary_to_master_table_ORFeome.R") # R subscript
}

# Call ORFtag hits -------------------------
if(F)
{
  file.edit("R/callORFtagHits.R")
  file.edit("R/callORFtagHitsStrandBias.R")
}

# Plots --------------------------------
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
  file.edit("R/vl_boxplot.R")

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

# Linear models evaluation -----------------
if(F)
{
  # Extract model performance stat
  file.edit("R/getModelEquation.R")
  file.edit("R/getModelRMSErsq.R")
  file.edit("R/getModelExpVar.R")
}

# STRING db tools --------------------------
if(F)
{
  # Get interactions
  file.edit("R/stringGetDB.R")
  file.edit("R/stringInteraction.R")
  file.edit("R/stringToIgraph.R")

  # Methods to plot
  file.edit("R/plot.vl_STRING.R")
}

# Miscellaneous ----------------------------
if(F)
{
  file.edit("R/gaussianBlur.R")
  file.edit("R/plot_enhancer_motifs.R") # Not reviewed
  file.edit("R/alignSanger.R") # Not reviewed but contains multiple alignment function
  file.edit("R/select_actMatched_controls.R") # Not reviewed, but contains code usable to sample closest observation
}
