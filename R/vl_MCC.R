#' Compute Matthew's Correlation Coefficient
#'
#' Matthew’s correlation coefficient (MCC) is a balanced binary-classification metric that measures
#' the correlation between predicted and true labels.
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#'
#' @return Matthew's correlation coefficient (ranges from -1 to 1, see examples).
#' 
#' @examples
#' # example code
#' vl_MCC(predicted= c(1,1,1,0,0), label= c(T,T,T,F,F))
#' vl_MCC(predicted= c(1,1,1,0,0), label= !c(T,T,T,F,F))
#' 
#' @export
vl_MCC <- function(predicted, label)
{
  # Check
  if(!is.logical(label))
    label <- as.logical(label)
  
  # Compute mcc
  yardstick::mcc_vec(
    truth    = factor(label, c(TRUE, FALSE)),
    estimate = factor(predicted > 0.5, c(TRUE, FALSE)),
    na_rm= TRUE
  )
}