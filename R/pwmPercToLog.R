#' PWM perc to log2
#'
#' @param perc_pwm A percentage PWM, where all columns sum to 1.
#' @param check.colsums If set to TRUE (default), ensures that all the columns of the perc_pwm sum to one.
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
                         check.colsums= TRUE,
                         pseudocount= 1e-5,
                         bg= c("A"= .25, "C"= .25, "G"= .25, "T"= .25))
{
  # Checks
  if(check.colsums && !all(colSums(perc_pwm)==1))
    stop("All the columns of perc_pwm should sum to 1.")
  # Make background matrix
  bg <- as.matrix(bg)
  bg <- lapply(seq(ncol(perc_pwm)), function(x) bg)
  bg <- do.call(cbind, bg)
  bg <- bg[c("A", "C", "G", "T"),]
  # Compute log2 odd ratio
  log2(perc_pwm/bg+pseudocount)
}
