#' importJASPAR
#'
#' Imports a .pfm text file containing JASPAR PFMs and parse them into lists of PPM and PWM TFBSTools::PWMatrix objects.
#'
#' @param combinedFile Path to a combined .pfm file containing JASPAR-formatted raw PFMs.
#' @param pseudocount Pseudocount to be added to PPMs before logging. Default= 1e-5.
#' @param bg Background probabilities. Default= c("A"= .25, "C"= .25, "G"= .25, "T"= .25).
#'
#' @return A list containing PPM and PWM matrices stored in TFBSTools::PWMatrix sublists.
#'
#' @examples
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' importJASPAR(pfm.file)
#'
#' @export
importJASPAR <- function(combinedFile,
                         pseudocount= 1e-5,
                         bg= c("A"= .25, "C"= .25, "G"= .25, "T"= .25))
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
  name <- gsub("^[^\\.]*\\.[^\\.]*\\.(.*)$", "\\1", IDs)
  if(anyDuplicated(name)){
    warning("Duplicated motif names were made unique.")
    name <- make.unique(name)
  }
  IDs <- gsub("^([^\\.]*\\.[^\\.]*)\\..*$", "\\1", IDs)
  if(anyDuplicated(IDs)){
    warning("Duplicated motif IDs were made unique.")
    IDs <- make.unique(IDs)
  }

  # Save PPM in PWMatrix object ----
  pwms_perc  <- lapply(seq(mat), function(i) {
    .c <- apply(mat[[i]], 2, function(x) x/sum(x))
    TFBSTools::PWMatrix(ID= IDs[i],
                        name = name[i],
                        profileMatrix = .c)
  })
  names(pwms_perc) <- name

  # log2 ratio PWMatrix object ----
  pwms_log_odds <- lapply(seq(pwms_perc), function(i) {
    .c <- pwms_perc[[i]]@profileMatrix
    .c <- pfmToPWM(matrix = .c,
                   pseudocount = pseudocount,
                   bg = bg)
    TFBSTools::PWMatrix(ID= IDs[i],
                        name = name[i],
                        profileMatrix = .c)
  })
  names(pwms_log_odds) <- name

  # Return perc pwm ---
  res <- list(pwms_perc= pwms_perc,
              pwms_log_odds= pwms_log_odds)
  return(res)
}

