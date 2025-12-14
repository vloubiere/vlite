setwd("/groups/stark/vloubiere/vlite/")
devtools::load_all("./")

# Tests/wip --------------------------------
if(F) {
  file.edit("inst/test/newUmiCollapsingStrat.R")
  # SOM clustering wrapper
  file.edit("R/somClustering.R")
  file.edit("R/somClusteringHelperFunctions.R")
  # PRO-seq pipeline with ncRNA decoy (test)
  file.edit("R/create_rRNA_tRNA_bowtie_index_mm10.R")
  file.edit("R/proseqProcessing_ncRNAdecoy.R")
}

# PWM tools --------------------------------
if(F) {
  # Switch between normalizations
  file.edit("R/pfmToPWM.R") # pfm to PWM
  file.edit("R/pwmToICM.R") # PWM to ICM

  # Import JASPAR
  file.edit("R/importJASPAR.R") # Import a combined JASPAR file

  # Plot DNA letters/logos
  file.edit("R/plotDNAletter.R") # Helper funciton
  file.edit("R/vl_seqLogo.R") # Plot the logo of a given PWM
  file.edit("R/vl_percLogo.R") # Plot percentage or frequency matrix
  file.edit("R/addSeqLogo.R") # Add PWM to an existing plot
  file.edit("R/addMotifs.R") # Add motifs to a heatmap
}

# AUC and model performance ----------------
if(F) {
  file.edit("R/vl_rocAUC.R") # Compute roc AUC and NES
  file.edit("R/vl_rocAUC_grid.R") # Parallelized version
  file.edit("R/vl_PPV.R") # Used in Shenzhi's paper
  file.edit("R/vl_PR_AUC.R") # Precision Recall AUC (imbalanced datasets)
  file.edit("R/vl_MCC.R") # Matthew's correlation coeff for binary labels
  file.edit("R/vl_TPR.R") # TRUE positive rate curve
}

# Single-cell tools ------------------------
if(F) {
  # Reads processing
  file.edit("R/cmd_alignCellRanger.R") # align sn-RNA-seq (GEX)
  file.edit("R/cmd_cellRangerAggr.R") # Merge multiple 10x runs (used to merge GEX)
  file.edit("R/cmd_alignCellRangerArc.R") # align multiome data (GEX+ATAC)

  # Velocity (GEX)
  file.edit("R/cmd_velocyto.R") # Run velocyto on GEX
  file.edit("inst/Rscript/filter_loom_file.R") # Function
  file.edit("R/cmd_filter_loom_file.R") # Wrapper

  # My SCENIC ----
  file.edit("inst/Rscript/infer_candidate_regulons.R") # Function
  file.edit("R/cmd_infer_candidate_regulons.R") # Wrapper (for single gene)
  file.edit("R/SCENIClite_infer_regulons.R") # All genes in parallel

  # sc-RNA-Seq
  file.edit("R/sc_computeMarkerGenes.R") # Identify marker genes from seurat object cluster(s)
  file.edit("R/sc_topMarkers.R") # Select top marker genes
  file.edit("R/sc_markerHeatmap.R") # Plot and cluster top marker genes from single-cell transcriptome data
  file.edit("R/sc_colors.R") # Default colors
  file.edit("R/sc_UMAP.R") # Plot rasterized UMAP from Seurat onbject
}

# bedtools ---------------------------------
if(F) {
  file.edit("R/importBed.R") # Import bed file
  file.edit("R/resizeBed.R") # Resize regions
  file.edit("R/binBed.R") # Bin each region
  file.edit("R/collapseBed.R") # Collapse overlapping regions
  file.edit("R/closestBed.R") # Find the closest regions (typically between two sets)
  file.edit("R/covBed.R") # Compute the number of overlaps between two sets of regions
  file.edit("R/overlapBed.R") # Compute overlaps between two sets if regions
  file.edit("R/intersectBed.R") # Compute intersection between two sets of regions
  file.edit("R/subtractBed.R") # Subtract a set of regions to another set of regions
  file.edit("R/clipBed.R") # Clip a set of regions using a second set of regions
  file.edit("R/exportBed.R") # Export bed file

  file.edit("R/randomRegionsBed.R") # Select a random set of regions from a (typically larger) set or region

  file.edit("R/enrichBed.R") # Assess whether the overlap between two set of regions is significanly large than expected
}

# Gene ontologies --------------------------
if(F) {
  file.edit("R/vl_GOenrich.R") # Compute GO enrichment for a set (or clusters) of Drosophila/mouse/human genes.
}

# bamtools ---------------------------------
if(F) {
  # Import in R
  file.edit("R/importBam.R") # Import bam file in R (Rsamtools wrapper)
  file.edit("R/covBam.R") # Wrapper around bedtools coverage
}

# BSgenomes tools -------------------------
if(F) {
  file.edit("R/getBSsequence.R") # Wrapper around ?getSeq
  file.edit("R/getBSgenomeSize.R") # Wrapper around ?seqinfo
  file.edit("R/randomRegionsBSgenome.R") # Uses ?randomRegionsBed
}

# fqtools ---------------------------------
if(F) {
  file.edit("inst/perl/parseFq.pl") # Function
  file.edit("R/importFq.R") # ?fread wrapper to quickly import fq files in R
}

# Download/upload data ---------------------
if(F) {
  # VBCF
  file.edit("inst/shell/download_bam_vbc.sh") # Download VBC files from fsk3
  # Dropbox API
  file.edit("R/dropboxUpload.R")
  file.edit("R/dropboxDownload.R")
  # SRA toolkit wrapper
  file.edit("R/cmd_downloadSRA.R") # Generates a command to download SRA files
  # fqs (or other...)
  file.edit("R/cmd_download.R") # A wget wrapper to download files
}

# SLURM wrappers ---------------------------
if(F) {
  file.edit("R/bsub.R") # bsub gridengine from Stark lab
  file.edit("R/vl_submit.R") # Wrapper around bsub to automatically submit pipeline commands
  file.edit("R/vl_squeue.R") # See user jobs
  file.edit("R/vl_scancel.R") # Cancel all jobs except interactive R session
  file.edit("R/vl_last_err.R") # Show the most recent error file from a logs folder
}

# Motifs analyses tools --------------------
if(F) {
  # Reverse complement DNA sequence
  file.edit("R/revCompDNA.R")

  # Motif counts
  file.edit("R/vl_motifCounts.R") # Wrapper around ?motifmatchr::matchMotifs to count motifs

  # Motif positions
  file.edit("R/vl_motifPos.R") # Wrapper around ?motifmatchr::matchMotifs to map motifs
  file.edit("R/motifPosToBed.R") # motiPos to bed
  file.edit("R/motifPosToMatrix.R") # motiPos to matrix (heatmap)

  # Motif enrichment (based on counts)
  file.edit("R/vl_motifEnrich.R") # Identify motifs over-represented in a set of sequences

  # Lasso regression motifs
  file.edit("R/motifLassoRegression.R") # Predict a reponse variable using motif counts

  # Download iCisTarget output
  file.edit("R/download_iCistarget.R") # Helper function to download iCistarget output
}

# bwtools ---------------------------------
if(F) {
  # Check content
  file.edit("R/bwGetSeqlengths.R") # Retrieve sequence lengths

  # Compute coverage
  file.edit("R/bwCoverage.R")
  file.edit("R/bwBinnedCoverage.R") # First bin input regions and quantifies (heatmap...)

  # Bigwig screenshot
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
  file.edit("R/importMetadataSheet.R") # Deprecated

  # Demultiplexing commands -------
  file.edit("inst/perl/vbc_tar_extract_head.pl") # Check the reads (for debugging)
  file.edit("R/checkVBCfile.R") # Wrapper
  file.edit("inst/perl/vbc_tar_demultiplexing.pl") # Perl subscript
  file.edit("inst/perl/vbc_bam_demultiplexing.pl") # Perl subscript
  file.edit("R/cmd_demultiplexVBCfile.R") # Wrapper

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
  # Pipeline ----
  file.edit("R/rnaseqProcessing.R")
  # Rsubread alignment (RNA-Seq)
  file.edit("R/cmd_alignRnaRsubread.R") # Wrapper
  file.edit("inst/Rscript/align_rna_Rsubread.R") # R subscript
  # Count reads (Rsubread)
  file.edit("R/cmd_countRsubread.R") # Wrapper
  file.edit("inst/Rscript/count_Rsubread.R") # R subscript
  # DESeq2 ----
  file.edit("R/cmd_DESeq2.R") # Wrapper
  file.edit("inst/Rscript/DESeq2_analysis.R") # R subscript

  # PRO-Seq -----------------------
  # Pipeline ----
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
  # Compute counts ----
  file.edit("inst/Rscript/create_PROseq_annotations.R") # Create the mm10 annotation files
  file.edit("inst/Rscript/create_PROseq_blacklisted_regions.R") # Create a mask containing tRNA coordinates...
  file.edit("inst/Rscript/count_PROseq_reads.R") # R subscript
  file.edit("R/cmd_countPROseqReads.R") # Wrapper
  # DESeq2 analysis ----
  file.edit("inst/Rscript/DESeq2_PROseq_analysis.R") # R subscript
  file.edit("R/cmd_DESeq2_PROseq.R") # Wrapper

  # STAP-Seq -----------------------
  # Pipeline (wrapper around the PRO-Seq) ----
  file.edit("R/stapseqProcessing.R")

  # CUT&RUN -----------------------
  # Pipeline ----
  file.edit("R/cutnrunProcessing.R")
  # QC
  file.edit("inst/Rscript/cutnrunQC.R") # R subscript
  file.edit("R/cmd_cutnrunQC.R") # Wrapper
  # MACS2 peak calling ----
  file.edit("R/cmd_peakCalling.R") # Wrapper
  # Confident peaks
  file.edit("R/cmd_confidentPeaks.R") # Wrapper
  file.edit("inst/Rscript/confident_peaks.R") # R subscript
  # Home made peak caller
  file.edit("R/cmd_peakCallingFromBw.R") # Wrapper
  file.edit("inst/Rscript/peakCallingFromBw.R") # R subscript

  # ORFtag -------------------------
  # Pipeline ----
  file.edit("R/orftagProcessing.R")
  file.edit("R/orftagQC.R")
  # Collapse bam unique insertions
  file.edit("R/cmd_collapseBam.R")
  # Assign reads to exons
  file.edit("inst/Rscript/create_ORFtag_annotations.R") # non-first exon gtf annotations
  file.edit("R/cmd_assignInsertions.R") # Wrapper
  file.edit("inst/Rscript/assign_ORFtag_insertions.R") # R subscript
  # Call hits ----
  file.edit("R/callOrfTagHits.R") # Wrapper
  file.edit("R/callOrftagHitsStrandBias.R") # R subscript

  # ORFeome -----------------------
  # Assemble dictionary ----
  file.edit("inst/Rscript/ORFeome_make_dictionary.R") # R subscript
  file.edit("R/cmd_orfeomeDictionary.R") # Wrapper
  # Pipeline ----
  file.edit("R/orfeomeProcessing.R")
  # Count ORFeomes barcode (ORFeome)
  file.edit("inst/Rscript/BC_counts.R") # R subscript
  file.edit("R/cmd_countBCreads.R") # Wrapper
  # Call hits MAGECK
  file.edit("inst/Rscript/compute_MAGECK_count_tables_ORFeome.R") # R subscript
  file.edit("inst/Rscript/merge_gene_summary_to_master_table_ORFeome.R") # R subscript
  file.edit("inst/Rscript/volcano_plots_MAgECK.R") # R subscript
  file.edit("R/cmd_MAGECK_ORFeome.R") # Wrapper

  # STARR-Seq ---------------------
  # Pipeline ----
  file.edit("R/starrseqProcessing.R")
  # UMI collapsing
  file.edit("R/cmd_umiCountsFromBam.R") # Wrapper
  file.edit("inst/Rscript/umiCountsFromBam.R") # R subscript
  # log2 Ratio bigwig
  file.edit("R/cmd_logRatioBigwig.R") # Wrapper
  file.edit("inst/Rscript/logRatioBigwig.R") # R subscript
}

# Ploting methods ----------------------
if(F) {
  # Plot enrichment analyses output
  file.edit("R/plot_vl_enr.R")
  file.edit("R/plot_vl_enr_clusters.R")
  file.edit("R/vl_plot_auc_enrichment.R")
}

# Plot legends -------------------------
if(F) {
  # Graphical parameters
  file.edit("R/vl_par.R")

  # Legends
  file.edit("R/tiltAxis.R")
  file.edit("R/vl_legend.R")
  file.edit("R/heatkey.R")
  file.edit("R/balloons_key.R")

  # Extra labels
  file.edit("R/addRsq.R")
  file.edit("R/addPcc.R")
  file.edit("R/addPval.R")
}

# Plots --------------------------------
if(F) {
  # Scatterplots
  file.edit("R/vl_plot.R")
  file.edit("R/rasterScatterplot.R")
  file.edit("R/addDensity.R")

  # Add repel labels (imperfect, but the best I can do)
  file.edit("R/helperFunctionRepelLabels.R")
  file.edit("R/addRepelLabels.R")

  # Heatmap
  file.edit("R/helperFunctions_heatmap.R")
  file.edit("R/vl_heatmap.R")

  # Other plots
  file.edit("R/vl_boxplot.R")
  file.edit("R/vl_pie.R")
  file.edit("R/alluvial.R")
  file.edit("R/balloons_plot.R")
  file.edit("R/MAplot.R")
  file.edit("R/plotTable.R")
  file.edit("R/upsetPlot.R")
}

# Deeplearning tools -----------------------
if(F) {
  # Deep learning contrib
  file.edit("R/importContrib.R")
  file.edit("R/contribSeqlets.R")
  file.edit("R/contribSeqLogo.R")
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
  file.edit("R/vl_cache_file.R")
  file.edit("R/gaussianBlur.R")
  file.edit("R/plot_enhancer_motifs.R") # Not reviewed
  file.edit("R/alignSanger.R") # Not reviewed but contains multiple alignment function
  file.edit("R/select_actMatched_controls.R") # Not reviewed, but contains code usable to sample closest observation
}
