#' LASSO regression using TF motifs
#'
#' @param response Variable to predict
#' @param counts Matrix of motif counts
#'
#' @return Return a LASSO model
#'
#' @examples
#' test <- motifLassoRegression(ATACSeq_FC, counts)
#'
#' @export
motifLassoRegression <- function(response,
                                 counts)
{
  # Setting alpha = 1 implements lasso regression
  lambdas <- 10^seq(2, -3, by = -.1)
  lasso_reg <- cv.glmnet(counts,
                         response,
                         alpha = 1,
                         lambda = lambdas,
                         standardize = TRUE,
                         nfolds = 5)
  # Best  lambda
  lambda_best <- lasso_reg$lambda.min
  # Modelling
  model <- glmnet(counts,
                  response,
                  alpha = 1,
                  lambda = lambda_best,
                  standardize = TRUE)

  # Predicted values
  # predict_test <- predict(model,
  #                         s = lambda_best,
  #                         newx = counts)
  # Use to compute adjusted Rsquare. Idea:
  # Number of observations
  # n <- length(response)
  # Number of non-zero coefficients (excluding the intercept)
  # p <- sum(coef(model) != 0) - 1  # Subtract 1 for the intercept
  # Residual sum of squares
  # SSres <- sum((response - predict_test)^2)
  # Total sum of squares
  # SStot <- sum((response - mean(response))^2)
  # Adjusted R-squared
  # R2_adj <- 1 - (SSres / (n - p - 1)) / (SStot / (n - 1))

  # Return
  return(model)
}
