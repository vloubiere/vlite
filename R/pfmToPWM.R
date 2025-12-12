#' Frequency matrix to PWM
#'
#' Matrices used to represent DNA binding motifs commonly come in four flavors:
#'
#' - PFM (Position Frequency Matrix): Raw counts of each nucleotide at each position.
#'
#' - PPM (Position Probability Matrix): Column-normalized version of the PFM where
#'   counts are converted to probabilities (each column sums to 1).
#'
#' - PWM/PSSM (Position Weight/Specific Scoring Matrix): Log-odds scores that compare
#'   the PPM to a background model. Typically: log2(observed_probability / background_probability).
#'   These matrices are used for scoring sequences and calling motif hits.
#'
#' - ICM (Information Content Matrix): Visualization-oriented transformation derived
#'   from the PPM using Shannon information. For DNA, the per-position information
#'   content is I = 2 - H, where H is the column entropy in bits; base heights in
#'   sequence logos are p_base * I. The maximum per-position information is 2 bits.
#'
#' @param matrix A PPM or a PFM that will be converted to PFM before computing the output PWM.
#' @param pseudocount Pseudocount to be added to the PFM columns before logging. Default= 1e-5.
#' @param bg Background probabilities. Default= c("A"= .25, "C"= .25, "G"= .25, "T"= .25).
#'
#' @return Returns a PWM.
#' @export
#'
#' @examples
#' PPM <- matrix(
#'   c(
#'     0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,
#'     0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,
#'     0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,
#'     0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22
#'   ), nrow= 4
#' )
#' PWM <- pfmToPWM(PPM)
#'
pfmToPWM <- function(matrix,
                     pseudocount= 1e-5,
                     bg= c("A"= .25, "C"= .25, "G"= .25, "T"= .25))
{
  # Checks
  if(!is.numeric(pseudocount))
    stop("pseudocount should be numeric.")
  if(!is.numeric(bg) | !identical(c("A", "C", "G", "T"), names(bg)))
    stop("bg should be a numeric vector  with names 'A', 'C', 'G', 'T'")
  if(!all(colSums(matrix)==1)) {
    warning("matrix contains columns that do not sum to 1 -> re-normalizing to PFM.")
    matrix <- apply(matrix, 2, function(x) x/sum(x))
  }
  if(is.null(colnames(matrix)))
    colnames(matrix) <- seq(ncol(matrix))


  # Make background matrix
  bg <- as.matrix(bg)
  bg <- lapply(seq(ncol(matrix)), function(x) bg)
  bg <- do.call(cbind, bg)
  bg <- bg[c("A", "C", "G", "T"),]

  # Compute log2 odd ratio
  pwm <- log2(matrix/bg+pseudocount)
  rownames(pwm) <- rownames(bg)

  # Return
  return(pwm)
}
