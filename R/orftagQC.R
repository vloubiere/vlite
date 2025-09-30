#' Perform QC of ORFtag screen batches
#'
#' @param sampleIDs A vector of sample IDs.
#' @param align.stats A vector of paths to alignment statistics files matching the sampleIDs.
#' @param counts.same.strand A vector of path of assigned insertions files matching the sampleIDs.
#' @param dc.cutoff The minimum number of duplicate counts for an insertion to be considered. Default= 0
#' @param output.prefix The output prefix for the output files.
#' @param output.folder Folder where output files will be saved. Default= "db/QC/ORFtag/".
#' @param pdf.width Width of the output pdf in inches. Default= 7.
#' @param pdf.height Height of the output pdf in inches. Default= 7.
#'
#' @return
#' @export
#'
#' @examples
orftagQC <- function(sampleIDs,
                     align.stats,
                     counts.same.strand,
                     dc.cutoff= 0,
                     output.prefix,
                     output.folder= "db/QC/ORFtag/",
                     pdf.width= 9,
                     pdf.height= 9)
{
  # Checks ----
  if(any(!file.exists(align.stats)))
    stop("Some alignment files do not exist or could not be found.")
  if(any(!file.exists(counts.same.strand)))
    stop("Some counts files do not exist or could not be found.")

  # Metadata ----
  meta <- data.table(sampleID= sampleIDs,
                     align.stats,
                     counts.same.strand)
  meta <- unique(meta)

  # For each cutoff ----
  output.files <- list(counts_table= c(),
                       pdf.qc= c())
  for(min.dc in dc.cutoff) {

    # Create output file names ----
    dir.create(output.folder, showWarnings = F, recursive = T)
    counts.table <- file.path(output.folder, paste0(output.prefix, "_", min.dc, "_counts.txt"))
    QC.pdf <- file.path(output.folder, paste0(output.prefix, "_", min.dc, "_QC.pdf"))

    # Compute stats ----
    res <- meta[, {
      # Alignment stats ----
      align <- readLines(align.stats)
      total_reads <- as.integer(unlist(tstrsplit(align[1], " ", keep=1)))
      perc <- unlist(tstrsplit(rev(align)[1], " ", keep=1))

      # Insertions stats ----
      coll <- fread(counts.same.strand)
      uniq <- nrow(coll)
      coll <- coll[dist<2e5]
      assigned <- nrow(coll)
      assigned.perc <- paste0(round(assigned/uniq*100, 1), "%")
      coll <- coll[ins_cov>=min.dc]
      dc.count <- nrow(coll)
      dc.count.perc <- paste0(round(dc.count/assigned*100, 1), "%")

      # Insertions IDs ----
      ins <- unique(coll[, .(seqnames, start, end, strand)])
      ins <- ins[, paste0(seqnames, ":", start, "-", end, ":", strand)]

      # Print table ----
      .c <- data.table("Aligned reads"= formatC(total_reads, big.mark = ",", format = "d"),
                       "Aligned %"= perc,
                       "Unique insertions"= formatC(uniq, big.mark = ","),
                       "Assigned insertions"= formatC(assigned, big.mark = ","),
                       "Assigned %"= assigned.perc,
                       "Ins. DC >="= paste0(formatC(dc.count, big.mark = ","), " (", dc.count.perc, ")"),
                       ins.IDs= list(ins))
      setnames(.c, "Ins. DC >=", paste0("Ins. DC >=", min.dc))
    }, .(sampleID, align.stats, counts.same.strand)]
    res$align.stats <- res$counts.same.strand <- NULL

    # Save statistics ----
    fwrite(res[, !"ins.IDs"],
           counts.table,
           col.names = T,
           row.names = F,
           sep= "\t",
           quote= F)

    # Compute overlaps ----
    perc <- lapply(res$ins.IDs, function(x)
      sapply(res$ins.IDs, function(y)
        round(sum(x %in% y)/length(x)*100)))
    perc <- do.call(rbind, perc)
    rownames(perc) <- res$sampleID
    colnames(perc) <- res$sampleID
    diag(perc) <- NA

    # Initiate plot ----
    pdf(QC.pdf, pdf.width, pdf.height)

    # Plot counts table ----
    plotTable(res[, !"ins.IDs"], cex = .6)

    # Plot overlaps matrix  ----
    vl_par(mai= rep(3.2, 4))
    vl_heatmap(perc,
               cluster.rows = F,
               cluster.cols = F,
               show.numbers = perc,
               breaks= seq(0, max(c(max(perc, na.rm = TRUE), 10)), length.out= 21),
               col= c("white", "red"),
               legend.title = "% Overlaps",
               numbers.cex = .4,
               main= rev(names(res))[2])
    dev.off()

    # Save output files ----
    output.files$counts_table[as.character(min.dc)] <- counts.table
    output.files$pdf.qc[as.character(min.dc)] <- QC.pdf
  }

  # Return output files
  return(output.files)
}
