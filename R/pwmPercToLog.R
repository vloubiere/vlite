#' PWM perc to log2
#'
#' @param perc_pwm A percentage PWM, where all columns sum to 1.
#' @param pseudocount Pseudocount to be added before logging. Default= 1e-5.
#' @param bg Background probabilities. Default= c("A"= .25, "C"= .25, "G"= .25, "T"= .25).
#'
#' @return Returns a log2 odd ratio pwm
#' @export
#'
#' @examples
#' load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
#' pwm <- vl_Dmel_motifs_DB_full$pwms_perc[[5]]@profileMatrix
#' pwmPercToLog(pwm)
pwmPercToLog <- function(perc_pwm,
                         pseudocount= 1e-5,
                         bg= c("A"= .25, "C"= .25, "G"= .25, "T"= .25))
{
  # Checks
  if(!is.numeric(pseudocount))
    stop("pseudocount should be numeric.")
  if(!is.numeric(bg) | !identical(c("A", "C", "G", "T"), names(bg)))
    stop("bg should be a numeric vector  with names 'A', 'C', 'G', 'T'")

  # Make background matrix
  bg <- as.matrix(bg)
  bg <- lapply(seq(ncol(perc_pwm)), function(x) bg)
  bg <- do.call(cbind, bg)
  bg <- bg[c("A", "C", "G", "T"),]
  # Compute log2 odd ratio
  log2(perc_pwm/bg+pseudocount)
}
