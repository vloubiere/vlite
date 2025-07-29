#' Count Motif Occurrences in Genomic Regions
#'
#' @description
#' Counts motif occurrences in a set of sequences or genomic regions using a PWMatrixList.
#'
#' @param sequences A named character vector of sequences to analyze. This argument takes precedence over the bed argument.
#' @param bed Genomic ranges in a format compatible with ?importBed, from which genomic sequences will be retrieved when sequences is set to NULL.
#' @param pwm_log_odds A PWMatrixList (in log2 odds ratio format) containing motifs to count.
#' @param genome Genome to use as background when bg = "genome" and/or to retrieve sequences (when bed is specified).
#' Default= NULL.
#' @param bg Background model for motif detection. Options are "genome", "subject" (inferred from input sequences) or "even" (0.25, 0.25, 0.25, 0.25). Default= "genome".
#' @param p.cutoff p-value cutoff for motif detection. Default= 5e-5.
#' @param cleanup.cache Logical. If set to TRUE, clears cached intermediate results. Default= FALSE.
#'
#' @return A matrix of motif counts.
#'
#' @examples
#' # Download Dev enhancer from pe-STARR-Seq paper
#' tmp <- tempfile(pattern = ".xlsx")
#' download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-52921-2/MediaObjects/41467_2024_52921_MOESM4_ESM.xlsx",
#'               destfile = tmp)
#'
#' # Retrieve enhancers
#' dev <- readxl::read_xlsx(tmp, sheet = 2, skip = 1)
#' dev <- as.data.table(dev)
#' dev <- dev[group=="dev" & detail %in% c("medium", "strong")]
#'
#' # Load motifs and select JASPAR
#' load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds]
#'
#' # Compute counts from regions
#' mot1 <- vl_motifCounts(bed= dev,
#'                        genome= "dm3",
#'                        pwm_log_odds= sel)
#'
#' # Compute counts from sequences
#' mot2 <- vl_motifCounts(sequences = dev$enhancer_sequence,
#'                        genome= "dm3",
#'                        pwm_log_odds= sel)
#'
#' # Both approaches should be identical
#' identical(mot1, mot2)
#'
#' @export
vl_motifCounts <- function(sequences,
                           bed,
                           pwm_log_odds,
                           genome= NULL,
                           bg= "genome",
                           p.cutoff= 5e-5,
                           cleanup.cache= FALSE)
{
  # Checks ----
  if(!missing(sequences) && !missing(bed))
     warning("sequences are provided -> input bed will not be used.")
  if(missing(sequences) && missing(bed))
    stop("sequences of bed regions should be specified.")
  if(is.null(genome) && (missing(sequences) | bg=="genome"))
    stop("genome is set to NULL.")
  if(!"PWMatrixList" %in% class(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList, pwm_log_odds)
  if(!is.numeric(p.cutoff) || p.cutoff>1)
    stop("p.cutoff should be a numeric value <= 1")

  # Get sequences ----
  if(missing(sequences)) {
    bed <- importBed(bed)
    sequences <- getBSsequence(bed, genome)
  }

  # Use a temp directory for caching motif counts ----
  params <- list(pwm_log_odds,
                 sequences,
                 genome= genome,
                 p.cutoff= p.cutoff,
                 bg= bg)
  file.key <- digest::digest(params)
  file.cache <- file.path(tempdir(), paste0(file.key, ".rds"))

  # Compute counts ----
  if(cleanup.cache | !file.exists(file.cache)) {
    res <- motifmatchr::matchMotifs(pwms = pwm_log_odds,
                                    subject = sequences,
                                    genome = genome,
                                    p.cutoff = p.cutoff,
                                    bg = bg,
                                    out = "scores")@assays@data[["motifCounts"]]
    saveRDS(res, file.cache)
  } else
    res <- readRDS(file.cache)

  # Format ----
  res <- as.matrix(res)
  res <- as.data.table(res)
  setnames(res,
           TFBSTools::name(pwm_log_odds))

  # Save ----
  return(res)
}
