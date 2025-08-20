#' A function to analyse differential marker genes between two sets of clusters
#'
#' Compute pairwise comparison between clusters.
#'
#' @param dat A Seurat object.
#' @param select.clusters A character vector of the cluster names of interest.
#' @param exclude.clusters Specified clusters will be excluded before any analysis.
#' @param pairwise A character vector of cluster names to be used as controls. Default= TRUE,
#' @param output.prefix Basename of output files, which will also be used as title for the plots.
#' @param logFC.cutoff logFC cutoff (wilcoxon test). Default= 0.25.
#' @param wilcox.padj.cutoff p.adjust cutoff for the wilcoxon test p.value.
#' @param log2OR.cutoff log2OR cutoff (fisher.test). Default= 2.
#' @param fisher.padj.cutoff p.adjust cutoff for the fisher test p.value.
#' @param pdf.output.folder Folder where pdf files will be saved. Default= "pdf".
#' @param GO.output.folder Folder where GO.rds enrichment files will be saved. Default= "db/GO".
#' @param open.pdf Should the pdf file be opened At the end of the function? Default= TRUE.
#' @param cleanup.cache Should all enrichments be computed again, even if temporary files already exist? Default= FALSE.
#' @param verbose For debugging purposes. Default= FALSE.
#' @return
#' @export
#'
#' @examples
sc_diffMarkers <- function(
    dat,
    select.clusters,
    exclude.clusters= NULL,
    output.prefix,
    pairwise= T,
    logFC.cutoff= .25,
    wilcox.padj.cutoff= 0.05,
    log2OR.cutoff= 2,
    fisher.padj.cutoff= 1e-10,
    pdf.output.folder= "pdf/",
    GO.output.folder= "db/GO/",
    open.pdf= TRUE,
    cleanup.cache= FALSE,
    verbose= FALSE
) {
  # Checks ----
  stopifnot(class(dat)=="Seurat")
  if(isFALSE(pairwise))
    stop("pairwise should be TRUE or a character vector of control cluster names to be used as controls.")
  dir.create(pdf.output.folder, showWarnings = F, recursive = T)
  dir.create(GO.output.folder, showWarnings = F, recursive = T)

  # pdf file ----
  pdf.file <- if(!is.null(output.prefix)) {
    file.path(pdf.output.folder, paste0(output.prefix, ".pdf"))
  } else {
    pdf.file <- gsub(".pdf$", "", list.files(pdf.output.folder, "Rplot.*.pdf"))
    pdf.file <- make.unique(c(pdf.file, "Rplot"), sep = "")
    file.path(pdf.output.folder, paste0(rev(pdf.file)[1], ".pdf"))
  }
  print(paste0("file.show(", shQuote(pdf.file), ")"))

  # GO output file ----
  GO.output.file <- if(!is.null(output.prefix)) {
    file.path(GO.output.folder, paste0(output.prefix, ".rds"))
  } else {
    GO.output.file <- gsub(".rds$", "", list.files(pdf.output.folder, "GOtable.*.rds"))
    GO.output.file <- make.unique(c(pdf.file, "GOtable"), sep = "")
    file.path(GO.output.folder, paste0(rev(GO.output.file)[1], ".rds"))
  }
  print(paste0("readRDS(", shQuote(GO.output.file), ")"))

  # Compute enriched markers
  top <- sc_computeMarkerGenes(
    dat = dat,
    select.clusters = select.clusters,
    exclude.clusters = exclude.clusters,
    pairwise = pairwise,
    add.fisher = T,
    cleanup.cache = cleanup.cache,
    verbose= verbose
  )
  top[, diff.wilcox:= fcase(
    padj<wilcox.padj.cutoff & logFC>logFC.cutoff, "up",
    padj<wilcox.padj.cutoff & logFC<(-logFC.cutoff), "down",
    default = "unaffected"
  )]
  top[, diff.fisher:= fcase(
    padj.fisher<fisher.padj.cutoff & log2OR>log2OR.cutoff, "up",
    padj.fisher<fisher.padj.cutoff & log2OR<(-log2OR.cutoff), "down",
    default = "unaffected"
  )]

  # GO enrich ----
  if(cleanup.cache || !file.exists(GO.output.file)) {
    GO <- list()
    for(test in c("diff.wilcox", "diff.fisher")) {
      # Gene list
      genes_list <- list(
        up= top$symbol[top[[test]]=="up"],
        down= top$symbol[top[[test]]=="down"]
      )
      # Enrichment
      GO[[test]] <- vl_GOenrich(
        genes_list,
        geneUniverse.IDs = top$symbol,
        species = "Dm",
        keyType = "SYMBOL",
        cleanup.cache = cleanup.cache
      )
    }
    # Save
    saveRDS(GO, GO.output.file)
  } else
    GO <- readRDS(GO.output.file)

  # Plots ----
  while(dev.cur() > 1)
    dev.off()
  pdf(pdf.file, width= 18, height = 8)
  vl_par(mfrow= c(2,4))
  for(test in c("diff.wilcox", "diff.fisher")) {
    # Coordinates
    if(test=="diff.wilcox") {
      x <- top$logFC
      y <- -log10(top$padj)
      xlab <- "logFC"
      ylab <- "Wilcox padj (-log10)"
    }
    if(test=="diff.fisher") {
      x <- top$log2OR
      y <- -log10(top$padj.fisher)
      xlab <- "log2OR"
      ylab <- "Fisher padj (-log10)"
    }
    sig <- top[[test]]!="unaffected"
    col <- adjustcolor(ifelse(sig, "red", "lightgrey"), .3)
    # Plot
    par(mai = c(.9, 1.08, .9, 1.08))
    rasterScatterplot(
      x,
      y,
      col= col,
      xlab= xlab,
      ylab= ylab,
      main= paste0(test, ": ", output.prefix)
    )
    text(x[sig],
         y[sig],
         top$symbol[sig],
         cex= .4)
    par(mai = c(.9, 2.5, .9, 1.5))
    GO[[test]][, {
      if(any(padj<0.05 & log2OR.cutoff>0)) {
        plot.vl_enr_cl(.SD, top.enrich = 10, cex= .4, main= annotation)
      } else {
        plot.new()
        text(.5, .5, "No enrich")
      }
      print(annotation)
    }, keyby= annotation]
  }
  dev.off()

  # Open pdf
  if(open.pdf)
    file.show(pdf.file)

  # Return top table
  invisible(top)
}
