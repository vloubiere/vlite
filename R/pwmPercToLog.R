#' PWM perc to log2
#'
#' @param perc_pwm A percentage PWM, where all columns sum to 1
#' @param pseudocount Pseudocount. Default= 1e-5
#' @param bg Background probabilities. Default= c("A"= .25, "C"= .25, "G"= .25, "T"= .25).
#'
#' @return Returns a log2 odd ratio pwm
#' @export
#'
#' @examples
#' vl_Dmel_motifs_DB_full$pwms_perc[[5]]@profileMatrix
#' pwmPercToLog(vl_Dmel_motifs_DB_full$pwms_perc[[5]]@profileMatrix)
pwmPercToLog <- function(perc_pwm,
                         pseudocount= 1e-5,
                         bg= c("A"= .25, "C"= .25, "G"= .25, "T"= .25))
{
  bg <- as.matrix(bg)
  bg <- lapply(seq(ncol(perc_pwm)), function(x) bg)
  bg <- do.call(cbind, bg)
  bg <- bg[c("A", "C", "G", "T"),]
  log2(perc_pwm/bg+pseudocount)
}
