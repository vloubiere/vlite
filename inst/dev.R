setwd("/groups/stark/vloubiere/vlite/")
devtools::load_all("./")

# Tests/wip --------------------------------
if(F) {
  file.edit("inst/test/newUmiCollapsingStrat.R")
}

# bedtools ---------------------------------
if(F) {
  file.edit("R/bedTools.R")
  file.edit("R/randomRegionsBed.R")
  file.edit("R/enrichBed.R")
}

# Gene ontologies --------------------------
if(F) {
  file.edit("R/vl_GOenrich.R")
}

# bamtools ---------------------------------
if(F) {
  # Import in R
  file.edit("R/importBamRsamtools.R")
  file.edit("R/covBam.R")
}

# BSgenomes tools -------------------------
if(F) {
  file.edit("R/BSgenomeTools.R")
}

# fqtools ---------------------------------
if(F) {
  file.edit("inst/perl/parseFq.pl") # Function
  file.edit("R/importFq.R") # Wrapper
}

# Download/upload data ---------------------
if(F) {
  # Dropbox
  file.edit("R/dropboxUpload.R")
  file.edit("R/dropboxDownload.R")
  # SRA toolkit
  file.edit("R/cmd_downloadSRA.R")
  # fqs (or other...)
  file.edit("R/cmd_download.R")
}

# SLURM wrappers ---------------------------
if(F) {
  file.edit("R/bsub.R") # bsub gridengine
  file.edit("R/vl_submit.R")
  file.edit("R/vl_squeue.R")
  file.edit("R/vl_scancel.R")
  file.edit("R/vl_last_err.R")
}

# Motifs analyses tools --------------------
if(F) {
  # Reverse complement DNA sequence
  file.edit("R/revCompDNA.R")

  # pwm conversion
  file.edit("R/importJASPAR.R")
  file.edit("R/pwmPercToLog.R")

  # Motif counts and positions
  file.edit("R/vl_motifCounts.R")
  file.edit("R/vl_motifPos.R")
  file.edit("R/motifPosToMatrix.R")
  file.edit("R/motifPosToBed.R")

  # Motif enrichment
  file.edit("R/vl_motifEnrich.R")

  # Lasso regression motifs
  file.edit("R/motifLassoRegression.R")

  # Download iCisTarget output
  file.edit("R/download_iCistarget.R")
}

# bwtools ---------------------------------
if(F) {
  # Check content
  file.edit("R/bwGetSeqlengths.R")

  # Coverage
  file.edit("R/bwCoverage.R")
  file.edit("R/bwBinnedCoverage.R")

  # Screenshot
  file.edit("R/helperFunctions_bwScreenshot.R")
  file.edit("R/bwScreenshot.R")

  # Average tracks and heatmap
  file.edit("R/bwAverageTrack.R")
  file.edit("R/plot.bwAverageTrack.R")
  file.edit("R/bwHeatmap.R")
}

# Genomics pipelines -----------------------
if(F) {
  # Helper functions --------------
  file.edit("R/importMetadataSheet.R") # NOT USED

  # Demultiplexing commands -------
  file.edit("R/cmd_demultiplexVBCfile.R") # Wrapper
  file.edit("inst/perl/vbc_tar_demultiplexing.pl") # Perl subscript
  file.edit("inst/perl/vbc_bam_demultiplexing.pl") # Perl subscript

  # General purpose commands -------
  # Trimming
  file.edit("R/cmd_trimIlluminaAdaptors.R") # Illumina adaptors
  # Alignment
  file.edit("R/cmd_alignBowtie.R") # Bowtie1 alignment
  file.edit("R/cmd_alignBowtie2.R") # Bowtie2 alignment
  # Merge several bigwig files
  file.edit("R/cmd_mergeBigwig.R") # Wrapper
  # Bedgraph to bigwig
  file.edit("R/cmd_bedgraphToBigwig.R") # Wrapper
  file.edit("inst/Rscript/bedgraph_to_bigwig.R") # R subscript
  # BAM to bigwig
  file.edit("R/cmd_bamToBigwig.R") # Wrapper
  file.edit("inst/Rscript/bam_to_bigwig.R") # R subscript
  # BED to bigwig
  file.edit("R/cmd_bedToBigwig.R") # Wrapper
  file.edit("inst/Rscript/bed_to_bigwig.R") # R subscript

  # RNA-Seq -----------------------
  # Pipeline
  file.edit("R/rnaseqProcessing.R")
  # Rsubread alignment (RNA-Seq)
  file.edit("R/cmd_alignRnaRsubread.R") # Wrapper
  file.edit("inst/Rscript/align_rna_Rsubread.R") # R subscript
  # Count reads (Rsubread)
  file.edit("R/cmd_countRsubread.R") # Wrapper
  file.edit("inst/Rscript/count_Rsubread.R") # R subscript
  # DESeq2
  file.edit("R/cmd_DESeq2.R") # Wrapper
  file.edit("inst/Rscript/DESeq2_analysis.R") # R subscript

  # PRO-Seq -----------------------
  # Pipeline
  file.edit("R/proseqProcessing.R")
  # Trimming
  file.edit("R/cmd_trimProseqAdaptors.R") # Custom adaptors (PROseq)
  # Extract unaligned reads from bam
  file.edit("R/cmd_extractUnalignedReadsFromBam.R")
  # UMI counts
  file.edit("R/cmd_umiCollapsingProseq.R") # Wrapper
  file.edit("inst/Rscript/umiCollapsingProseq.R") # R subscript
  # PROSeq UMI counts to bigwig
  file.edit("R/cmd_umiToBigwigProseq.R") # Wrapper
  file.edit("inst/Rscript/umiToBigwigProseq.R") # R subscript
  # Count reads from .rds annotation file
  file.edit("R/cmd_countPROseqReads.R") # Wrapper
  file.edit("inst/Rscript/count_PROseq_reads.R") # R subscript
  # DESeq2
  file.edit("R/cmd_DESeq2_PROseq.R") # Wrapper
  file.edit("inst/Rscript/DESeq2_PROseq_analysis.R") # R subscript

  # CUT&RUN -----------------------
  # Pipeline
  file.edit("R/cutnrunProcessing.R")
  # QC
  file.edit("inst/Rscript/cutnrunQC.R") # R subscript
  file.edit("R/cmd_cutnrunQC.R") # Wrapper
  # MACS2 peak calling
  file.edit("R/cmd_peakCalling.R") # Wrapper
  # Confident peaks
  file.edit("R/cmd_confidentPeaks.R") # Wrapper
  file.edit("inst/Rscript/confident_peaks.R") # R subscript
  # Home made peak caller
  file.edit("R/cmd_peakCallingFromBw.R") # Wrapper
  file.edit("inst/Rscript/peakCallingFromBw.R") # R subscript

  # ORFtag -------------------------
  # Pipeline
  file.edit("R/orftagProcessing.R")
  file.edit("R/orftagQC.R")
  # Collapse bam unique insertions
  file.edit("R/cmd_collapseBam.R")
  # Assign reads to exons
  file.edit("R/cmd_assignInsertions.R") # Wrapper
  file.edit("inst/Rscript/assign_ORFtag_insertions.R") # R subscript
  # Call hits
  file.edit("R/callOrfTagHits.R") # Wrapper
  file.edit("R/callOrftagHitsStrandBias.R") # R subscript

  # ORFeome -----------------------
  # Assemble dictionary
  file.edit("R/orfeomeDictionary.R")
  # Pipeline
  file.edit("R/orfeomeProcessing.R")
  # Count ORFeomes barcode (ORFeome)
  file.edit("R/cmd_countBCreads.R") # Wrapper
  file.edit("inst/Rscript/BC_counts.R") # R subscript
  # Call hits MAGECK
  file.edit("R/cmd_MAGECK_ORFeome.R") # Wrapper
  file.edit("inst/Rscript/compute_MAGECK_count_tables_ORFeome.R") # R subscript
  file.edit("inst/Rscript/merge_gene_summary_to_master_table_ORFeome.R") # R subscript
  file.edit("inst/Rscript/volcano_plots_MAgECK.R") # R subscript

  # STARR-Seq ---------------------
  # Piepline
  file.edit("R/starrseqProcessing.R")
  # UMI collapsing
  file.edit("R/cmd_umiCountsFromBam.R") # Wrapper
  file.edit("inst/Rscript/umiCountsFromBam.R") # R subscript
  # log2 Ratio bigwig
  file.edit("R/cmd_logRatioBigwig.R") # Wrapper
  file.edit("inst/Rscript/logRatioBigwig.R") # R subscript

  # sc-RNA-Seq ---------------------
  file.edit("R/cmd_alignCellRanger.R")
}

# Plots --------------------------------
if(F) {
  # Plotting parameters
  file.edit("R/vl_par.R")
  file.edit("R/vl_plot.R")

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
  file.edit("R/helperFunctions_heatmap.R")
  file.edit("R/vl_heatmap.R")

  # Alluvial plot
  file.edit("R/alluvial.R")

  # Balloons plot
  file.edit("R/balloons_key.R")
  file.edit("R/balloons_plot.R")

  # Plot DNA letters/logos
  file.edit("R/plotDNAletter.R") # Used to plot logos
  file.edit("R/addSeqLogo.R")
  file.edit("R/vl_seqLogo.R")
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
  file.edit("R/addDensity.R")

  # Add repel labels (imperfect, but the best I can do)
  file.edit("R/helperFunctionRepelLabels.R")
  file.edit("R/addRepelLabels.R")

  # Upset plot
  file.edit("R/upsetPlot.R")
}

# Deeplearning tools -----------------------
if(F) {
  # Deep learning contrib
  file.edit("R/importContrib.R")
  file.edit("R/contribSeqlets.R")
  file.edit("R/contribSeqLogo.R")

  # Modelling diagnostics
  file.edit("R/vl_ROC_AUC.R")
  file.edit("R/vl_PR_AUC.R")
  file.edit("R/vl_PPV.R")
  file.edit("R/vl_mPCC.R")
  file.edit("R/vl_TPR.R")
}

# Linear models evaluation -----------------
if(F) {
  # Extract model performance stat
  file.edit("R/getModelEquation.R")
  file.edit("R/getModelRMSErsq.R")
  file.edit("R/getModelExpVar.R")
}

# STRING db tools --------------------------
if(F) {
  # Get interactions
  file.edit("R/stringGetDB.R")
  file.edit("R/stringInteraction.R")
  file.edit("R/stringToIgraph.R")

  # Methods to plot
  file.edit("R/plot.vl_STRING.R")
}

# Miscellaneous ----------------------------
if(F) {
  file.edit("R/gaussianBlur.R")
  file.edit("R/plot_enhancer_motifs.R") # Not reviewed
  file.edit("R/alignSanger.R") # Not reviewed but contains multiple alignment function
  file.edit("R/select_actMatched_controls.R") # Not reviewed, but contains code usable to sample closest observation
}
