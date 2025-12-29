#' importJASPAR
#'
#' Imports a .pfm text file containing JASPAR PFMs, and uses functions from the TFBStools package
#' to compute the following raw and normalized position matrices:
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
#' @param combinedFile Path to a combined .pfm file containing JASPAR-formatted raw PFMs.
#' @param bg Background probabilities. Default= c("A"= .25, "C"= .25, "G"= .25, "T"= .25).
#' @param pseudocount Pseudocount added to PFMs before computing PPMs, PWMs and ICMs. Default= .8.
#' @param schneider Should schneider correction be used for ICMs? Default= TRUE.
#'
#' @return A list containing PPM and PWM matrices stored in TFBSTools::PWMatrix sublists.
#'
#' @examples
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' importJASPAR(pfm.file)
#'
#' @export
importJASPAR <- function(combinedFile,
                         bg= c("A"= .25, "C"= .25, "G"= .25, "T"= .25),
                         pseudocount= 0.8,
                         schneider= TRUE)
{
  # Import JASPAR format ----
  lines <- readLines(combinedFile)

  # Import JASPAR format ----
  idx <- grep("^>", lines)
  mat <- lapply(idx, function(i) {
    mat <- do.call(
      rbind,
      lapply(1:4, function(j) {
        x <- unlist(tstrsplit(lines[i+j], " "))
        x <- x[x!=""]
        x <- as.numeric(x[3:(length(x)-1)])
      })
    )
    rownames(mat) <- c("A", "C", "G", "T")
    colnames(mat) <- NULL
    return(mat)
  })

  # IDs and names ----
  IDs <- unlist(tstrsplit(lines[idx], "\t", keep= 2))
  mot.names <- gsub("^[^\\.]*\\.[^\\.]*\\.(.*)$", "\\1", IDs)
  if(anyDuplicated(mot.names)){
    warning("Duplicated motif names were made unique.")
    mot.names <- make.unique(mot.names)
  }
  IDs <- gsub("^([^\\.]*\\.[^\\.]*)\\..*$", "\\1", IDs)
  if(anyDuplicated(IDs)){
    warning("Duplicated motif IDs were made unique.")
    IDs <- make.unique(IDs)
  }

  # Make PFM object ----
  PFM  <- lapply(seq(mat), function(i) {
    TFBSTools::PFMatrix(ID= IDs[i],
                        name = mot.names[i],
                        profileMatrix = mat[[i]],
                        bg = bg)
  })
  PFM <- do.call(TFBSTools::PFMatrixList, PFM)

  # Convert to PPM ----
  PPM  <- TFBSTools::toPWM(x= PFM,
                           type = "prob",
                           pseudocounts = pseudocount,
                           bg = bg)

  # Convert to PWM ----
  PWM  <- TFBSTools::toPWM(x= PFM,
                           type = "log2probratio",
                           pseudocounts = pseudocount,
                           bg = bg)

  # Convert to ICM ----
  ICM  <- TFBSTools::toICM(x= PFM,
                           pseudocounts = pseudocount,
                           schneider = schneider,
                           bg = bg)

  # Return perc pwm ---
  res <- list(
    metadata= data.table(motif= IDs, name= mot.names),
    PFM= PFM,
    PPM= PPM,
    PWM= PWM,
    ICM= ICM
  )
  return(res)
}

