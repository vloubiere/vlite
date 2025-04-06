#' Compute RMSE and R2 from observed and predicted values
#'
#' Comes from https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
#'
#' @param observed vector of observed values
#' @param predicted vector of predicted values
#'
#' @return RMSE and R2
#' @export
getModelRMSErsq <- function(observed, predicted)
{
  if(length(predicted)!=length(observed))
    print("Observed and predicted should have the same length")
  SSE <- sum((predicted - observed)^2)
  SST <- sum((observed - mean(observed))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/length(observed))

  # Model performance metrics
  data.table(RMSE = RMSE,
             Rsquare = R_square)
}
