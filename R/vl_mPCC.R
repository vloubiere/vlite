#' Compute Matthew's PCC
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#'
#' @return Mathhew's correlation coef
#' @export
vl_mPCC <- function(predicted,
                    label)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))

  # Make table
  conf_matrix <- table(pred= factor(predicted>0.5, c(FALSE, TRUE)),
                       obs= factor(as.logical(label), c(FALSE, TRUE)))
  TP <- as.numeric(conf_matrix[2, 2])
  TN <- as.numeric(conf_matrix[1, 1])
  FP <- as.numeric(conf_matrix[2, 1])
  FN <- as.numeric(conf_matrix[1, 2])
  # Compute
  mcPCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  # Return
  return(mcPCC)
}
