#' Motif enrichment analysis
#'
#' Compute motif enrichment between a set of regions and control regions
#'
#' @param counts A data.table or a list of data.tables containing motif counts for the regions of interest.
#' @param control.counts A data.table containing counts for control regions.
#' Should have the same number of columns as 'counts' data.table(s).
#' @param names A character vector or factor names matching control.counts' columns (i.e., motif IDs).
#' By default, column names are used as is.
#' @param fisher.alternative The alternative for the fisher test. By default, if counts is a data.table or
#' a list of length one, set to "two.sided". If counts contain at least two groups, set to "greater".
#' @param log2OR.pseudocount Numeric. A pseudocount added to the contingency table containing 0s, in order
#' to avoid infinite log2 odd ratios. Default= 0.5.
#' @param countFUN The function to be used for each motif, in order to compute the number of positive sequences.
#' By default, a sequence will be considered positive if it contains at least one motif.
#' Default= function(motifCount) sum(motifCount>=1).
#'
#' @examples
#' # Download Dev/Hk enhancers from pe-STARR-Seq paper
#' tmp <- tempfile(fileext = ".xlsx")
#' download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-52921-2/MediaObjects/41467_2024_52921_MOESM4_ESM.xlsx",
#'               destfile = tmp)
#'
#' # Retrieve enhancers
#' enh <- readxl::read_xlsx(tmp, sheet = 2, skip = 1)
#' enh <- as.data.table(enh)
#' enh <- enh[group=="dev" & detail %in% c("medium", "strong") | group=="hk"]
#'
#' # Negative controls
#' set.seed(1)
#' ctl <- controlRegionsBSgenome(bed = enh,
#'                               genome = "dm3",
#'                               no.overlap = TRUE)
#'
#' # Select JASPAR motifs
#' load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds]
#' pwms <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_perc]
#' mot.cluster <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_cluster]
#'
#' # Compute enrichment at developmental enhancers
#' counts <- vl_motifCounts(bed= enh[group=="dev"],
#'                          genome= "dm3",
#'                          pwm_log_odds= sel)
#' ctl.counts <- vl_motifCounts(bed= ctl,
#'                              genome= "dm3",
#'                              pwm_log_odds= sel)
#'
#' # Enrichment
#' enr <- vl_motifEnrich(counts,
#'                       control.counts = ctl.counts,
#'                       names = mot.cluster)
#'
#' # Collapse per motif cluster
#' coll <- enr[motif %in% enr[, motif[which.min(padj)], name]$V1]
#'
#' # Plot
#' vl_par(mai= c(.9, 3, .9, 1.2))
#' pl <- plot(coll)
#' addMotifs(plot.DT = pl,
#'           pwms = pwms)
#'
#' # Compare dev and Hk enhancers
#' counts1 <- vl_motifCounts(bed= enh,
#'                           genome= "dm3",
#'                           pwm_log_odds= sel)
#'
#' # Enrichment per cluster
#' enr1 <- vl_motifEnrich(split(counts1, enh$group),
#'                        control.counts = ctl.counts,
#'                        names = mot.cluster)
#'
#' # Collapse per motif cluster and plot
#' coll1 <- enr1[motif %in% enr1[, motif[which.min(padj)], name]$V1]
#'
#' # Plot
#' vl_par(mai= c(.9, 3, .9, 1.2))
#' pl <- plot(coll1)
#' addMotifs(plot.DT = pl,
#'           pwms= pwms)
#'
#' @return DT of enrichment values which can be plot using ?plot.vl_GO_enr
#' @export
vl_motifEnrich <- function(counts,
                           control.counts,
                           names= NULL,
                           fisher.alternative,
                           log2OR.pseudocount= .5,
                           countFUN= function(motifCount) sum(motifCount>=1))
{
  if(is.data.table(counts))
    counts <- list(set= counts)
  if(is.null(names(counts)))
    names(counts) <- seq(counts)
  if(anyDuplicated(names(counts)))
    stop("names(counts) should be unique!")
  if(any(!sapply(counts, is.data.table)))
    counts <- lapply(counts, as.data.table)
  if(!is.data.table(control.counts))
    control.counts <- as.data.table(control.counts)
  if(any(!sapply(counts, function(x) identical(colnames(x), colnames(control.counts)))))
    stop("counts / control.counts data.tables should contain the same column names")
  if(is.null(names))
    names <- colnames(control.counts)
  if(!is.character(names) && !is.factor(names))
    stop("names should be a vector of characters or of factors.")
  if(missing(fisher.alternative))
    fisher.alternative <- ifelse(length(counts)>1, "greater", "two.sided")

  # Melt control ----
  ctl <- melt(control.counts,
              measure.vars = names(control.counts),
              variable.name= "motif")
  ctl <- ctl[, .(ctl_hit= countFUN(value), ctl_total= .N), motif]

  # Add names ----
  add.names <- data.table(motif= colnames(control.counts),
                          name= names)
  ctl <- merge(ctl, add.names)

  # Melt counts ----
  cl <- lapply(counts, function(x) melt(x, measure.vars= names(x), variable.name= "motif"))
  cl <- rbindlist(cl, idcol = "cl")
  cl <- cl[, .(set_hit= sum(value>=1), set_total= .N), .(cl, motif)]

  # Compute enrichment ----
  enr <- merge(cl, ctl)
  enr[, c("OR", "pval"):= {
    # Contingency matrix
    mat <- c(set_hit, set_total - set_hit,
             ctl_hit, ctl_total - ctl_hit)
    mat <- matrix(mat, byrow= T, ncol= 2)
    # p.value
    .f <- fisher.test(mat,
                      alternative = fisher.alternative)
    # log2OR (pseudocount avoid Inf)
    if(any(mat==0)) {
      mat <- mat+log2OR.pseudocount
      .f$estimate <- (mat[1,1] * mat[2,2]) / (mat[2,1] * mat[1,2])
    }
    .(.f$estimate, .f$p.value)
  }, .(set_hit, set_total, ctl_hit, ctl_total)]

  # Compute log2OR and padj ----
  enr[, log2OR:= log2(OR)]
  enr[, padj:= p.adjust(pval, method = "fdr"), cl]
  enr$OR <- enr$pval <- NULL

  # Define class (for plotting methods) ----
  if(length(unique(enr$cl))>1) {
    setattr(enr, "class", c("vl_enr_cl", "data.table", "data.frame"))
  } else {
    setattr(enr, "class", c("vl_enr", "data.table", "data.frame"))
  }

  # Order and return enrichment table ----
  setcolorder(enr,
              c("cl", "name", "motif"))
  setorderv(enr,
            c("cl", "padj"))
  return(enr)
}
