#' reverseComplement
#'
#' Reverse complement DNA sequence character string.
#'
#' @param DNA_char Character vector, DNAString or DNAStringSet.
#' @param reverse Should the DNA sequence be reversed? Default= TRUE
#' @param complement Should the DNA sequence be complemented? Default= TRUE
#' @param as.character Should the sequences be returned as a character vector? Otherwise,
#' it will be returned as a DNAString or a DNAStringSet. Default= TRUE
#'
#' @examples
#' revCompDNA("ATCG")
#'
#' @return Reverse Complemented DNA sequence
#' @export
revCompDNA <- function(DNA_char,
                       reverse= TRUE,
                       complement= TRUE,
                       as.character= TRUE)
{
  # Checks
  if(is.character(DNA_char)) {
    DNA_char <- if(length(DNA_char==1)) {
      DNAString(DNA_char)
    } else
      DNAStringSet(DNA_char)
  }
  if(!class(DNA_char) %in% c("DNAString", "DNAStringSet"))
    stop("DNA_char should be a character vector, a DNAString or a DNAStringSet")

  # Reverse
  if(reverse)
    DNA_char <- IRanges::reverse(DNA_char)

  # Complement
  if(complement)
    DNA_char <- Biostrings::complement(DNA_char)

  # As character
  if(as.character)
    DNA_char <- as.character(DNA_char)

  # Return
  return(DNA_char)
}
