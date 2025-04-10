#' reverseComplement
#'
#' Reverse complement DNA sequence character string.
#'
#' @param DNA_char Input DNA character string.
#' @param complement Should the DNA sequence be complemented? Default= TRUE
#' @param reverse Should the DNA sequence be reversed? Default= TRUE
#'
#' @examples
#' revCompDNA("ATCG")
#'
#' @return Reverse Complemented DNA sequence
#' @export
revCompDNA <- function(DNA_char,
                       complement= TRUE,
                       reverse= TRUE)
{
  # Checks
  if(length(DNA_char)>1)
    stop("length DNA_char should be 1!")

  # Split letters
  .c <- strsplit(DNA_char, "")[[1]]

  # Check letters
  if(!all(.c %in% c("A", "T", "C", "G", "a", "t", "c", "g")))
    stop("All characters should be one of A, T, C, G, a, t, c, g")

  res <- sapply(.c, function(x)
  {
    switch(x,
           "A"= "T",
           "T"= "A",
           "C"= "G",
           "G"= "C",
           "a"= "t",
           "t"= "a",
           "c"= "g",
           "g"= "c")
  })
  if(reverse)
    res <- rev(res)
  return(paste0(res, collapse = ""))
}
