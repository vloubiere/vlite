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
#' @param combinedFile Path to a combined file containing JASPAR-formatted raw PFMs.
#' @param bg Background probabilities. Default= c("A"= .25, "C"= .25, "G"= .25, "T"= .25).
#' @param pseudocount Pseudocount added to PFMs before computing PPMs, PWMs and ICMs. Default= .8.
#' @param schneider Should schneider correction be used for ICMs? Default= FALSE (similar to
#' what is shown on JASPAR website).
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
                         schneider= FALSE,
                         simplify.names= T)
{
  # Import JASPAR combined file ----
  PFM  <- TFBSTools::readJASPARMatrix(combinedFile, matrixClass = "PFM")
  check.int <- sapply(PFM, function(x) is.integer(x@profileMatrix))
  if(!all(check.int))
    stop("All PFM matrices should contain integers.")
  
  # Simplify names ----
  if(simplify.names) {
    # Trim ID
    simp.names <- sapply(PFM, function(x) {
      gsub(paste0(x@ID, "."), "", x@name)
    })
    # Make unique
    if(anyDuplicated(simp.names))
      warning("Duplicated simplified names will be made unique")
    simp.names <- make.unique(simp.names)
    # Replace
    PFM <- lapply(seq_along(PFM), function(i) {
      x <- PFM[[i]]
      x@name <- simp.names[i]
      x
    })
    PFM  <- do.call(TFBSTools::PFMatrixList, PFM)
    names(PFM) <- TFBSTools::name(PFM)
  }

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

  # Final object ---
  res <- list(
    metadata= data.table(motif= TFBSTools::ID(PFM), name= TFBSTools::name(PFM)),
    PFM= PFM,
    PPM= PPM,
    PWM= PWM,
    ICM= ICM
  )
  
  # Return ---
  return(res)
}

