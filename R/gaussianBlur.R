#' Apply a Gaussian blur to a numeric vector
#'
#' Smooths a numeric vector using a Gaussian kernel.
#'
#' @param x Numeric vector containing signal to be blurred.
#' @param size Integer specifying the size of the Gaussian kernel (must be odd, default is 5).
#' @param sigma Numeric value specifying the standard deviation of the Gaussian kernel (default is 1).
#' @param padding_value Numeric calue used to pad the signal at both ends and avoid NAs (default is 0).
#'
#' @return A numeric vector of the same length as x, containing the blurred signal.
#' @export
gaussianBlur <- function(x,
                         size = 5,
                         sigma = 1,
                         padding_value = 0)
{
  # Checks
  if (size %% 2 == 0) stop("Size must be odd.")

  # Create Gaussian kernel
  kernel_x <- seq(-floor(size/2), floor(size/2), length.out = size)
  kernel <- exp(-0.5 * (kernel_x / sigma)^2)
  kernel <- kernel / sum(kernel)

  # Pad the vector
  Npads <- floor(size / 2)
  padded_x <- c(rep(padding_value, Npads), x, rep(padding_value, Npads))

  # Apply convolution
  blurred <- stats::filter(padded_x, kernel, sides = 2)

  # Remove padding
  blurred <- blurred[(Npads + 1):(length(blurred) - Npads)]
  return(as.numeric(blurred))
}
