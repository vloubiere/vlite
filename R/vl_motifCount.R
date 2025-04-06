#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start' and 'end' columns.
#' @param pwms_log_odds A PWMatrixList (in log2 odd ratio format) for which motif matches should be counted. Overrides sel and motifDB arguments (see above).
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome".
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff p.value cutoff used for motif detection. For enrichment analyses based on presence/absence of a motif, high cutoff might perform better (1e-4 or 5e-5) while for regression analyses, lower cutoffs might be preferred (5e-4). Default= 5e-5 (stringent).
#'
#' @examples
#' # Resize example peaks
#' SUHW <- vl_resizeBed(vl_SUHW_top_peaks, genome = "dm3")
#' STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, genome = "dm3")
#'
#' # Generate same number of random regions
#' random <- vl_control_regions_BSgenome(bed= STARR, genome= "dm3")
#'
#' # Count JAPSPAR motifs (see below to use custom list of PWMs)
#' suhw <- vl_motifCounts(SUHW, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' starr <- vl_motifCounts(top_STARR, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' ctl <- vl_motifCounts(random, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#'
#' # Starting from sequence instead of bed file
#' seq <- vl_getSequence(SUHW, genome= "dm3")
#' seq_suhw <- vl_motifCounts(seq, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' identical(suhw, seq_suhw)
#'
#' # Motifs can also be counted using a custom PWMatrixList, for example for promoter motifs:
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom))
#'   prom[[i]]@profileMatrix <- vl_pwm_perc_to_log2(prom_db$Pwms_perc[[i]]@profileMatrix)
#'
#' prom_motifs <- vl_motifCounts(STARR,
#'                                pwm_log_odds= prom,
#'                                genome= "dm3")
#'
#' @return Matrix of motif counts
#' @export
vl_motifCounts <- function(sequences, ...) UseMethod("vl_motifCounts")

#' @describeIn vl_motifCounts Method to extract sequences from BSgenome
#' @export
vl_motifCounts.data.table <- function(bed,
                                      genome,
                                      ...)
{
  sequences <- getBSsequence(bed, genome)
  vl_motifCounts.default(sequences,
                         genome= genome,
                         ...)
}

#' @describeIn vl_motifCounts Identify motifs in sequences
#' @export
vl_motifCounts.default <- function(sequences= NULL,
                                   pwm_log_odds,
                                   genome,
                                   bg= "genome",
                                   p.cutoff= 5e-5)
{
  # Checks ----
  if(!"PWMatrixList" %in% class(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList, pwm_log_odds)

  # Compute counts ----
  res <- motifmatchr::matchMotifs(pwm_log_odds,
                                  sequences,
                                  genome= genome,
                                  p.cutoff= p.cutoff,
                                  bg= bg,
                                  out= "scores")@assays@data[["motifCounts"]]
  res <- as.matrix(res)
  res <- as.data.table(res)
  setnames(res,
           TFBSTools::name(pwm_log_odds))

  # Save ----
  return(res)
}
