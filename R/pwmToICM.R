#' PWM to ICM (Information Content Matrix)
#'
#' Converts a DNA PWM (log2-odds vs background) into an Information Content Matrix (ICM).
#' Steps:
#'   1) Convert PWM to PPM using p_ij ∝ bg_i * 2^{PWM_ij}, normalized per column.
#'   2) Compute per-position information I_j = 2 - H_j (bits), where
#'      H_j = -sum_i p_ij * log2(p_ij).
#'   3) Return ICM with entries p_ij * I_j.
#'
#' @param pwm A numeric 4 x L matrix of log2-odds scores (PWM), rows named "A","C","G","T".
#' @param bg Background probabilities named c("A","C","G","T"). Default = uniform 0.25.
#'
#' @return A 4 x L ICM matrix (rows A,C,G,T), in bits.
#' @export
#'
#' @examples
#' # Example: start from a PPM and create a PWM, then convert to ICM
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
#' # convert to ICM
#' ICM <- pwmToICM(PWM)
pwmToICM <- function(PWM,
                     bg = c("A"=.25, "C"=.25, "G"=.25, "T"=.25))
{
  # Checks
  if (!is.matrix(PWM) || !is.numeric(PWM))
    stop("PWM must be a numeric matrix of log2-odds scores.")
  if(is.null(colnames(PWM)))
    colnames(PWM) <- seq(ncol(PWM))
  if (is.null(rownames(PWM)) ||
      !identical(rownames(PWM), c("A","C","G","T")))
    stop("Row names must be exactly c('A','C','G','T') in that order.")
  if (!is.numeric(bg) || !identical(names(bg), c("A","C","G","T")))
    stop("'bg' must be a numeric named vector with names c('A','C','G','T').")
  if (any(bg <= 0) || sum(bg)!=1)
    stop("'bg' must be positive and sum to 1.")

  # Convert PWM (log2 odds) to unnormalized log2 probabilities:
  # log2 p_ij ∝ log2 bg_i + PWM_ij
  log2p_unnorm <- sweep(PWM, 1, log2(bg), FUN = "+")

  # Normalize per column using base-2 softmax for numerical stability:
  # p_ij = 2^{log2p_unnorm_ij - log2sumexp_col}
  # Implement via natural exp: 2^x = exp(x * ln 2)
  ln2 <- log(2)
  shift <- apply(log2p_unnorm, 2, max)              # for stability
  M_shift <- sweep(log2p_unnorm, 2, shift, "-")
  E <- exp(M_shift * ln2)
  PPM <- sweep(E, 2, colSums(E), "/")

  # Guard against exact zeros for entropy
  PPM <- pmax(PPM, .Machine$double.eps)
  PPM <- sweep(PPM, 2, colSums(PPM), "/")  # re-normalize after clamping

  # Entropy per column and information (DNA: max = log2(4) = 2 bits)
  H <- -colSums(PPM * log2(PPM))
  I <- pmax(0, 2 - H)

  # Information Content Matrix
  ICM <- sweep(PPM, 2, I, "*")
  rownames(ICM) <- c("A","C","G","T")
  colnames(ICM) <- colnames(PWM)

  # Return
  return(ICM)
}
