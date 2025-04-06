#' Compute explained variance
#'
#' From Bernardo, computes the % explained variance for each predictor
#'
#' @param model An object containing the results returned by a model fitting function (e.g., lm or glm).
#' @return % of explained variance
#' @export
getModelexpVar <- function(model)
{
  af <- stats::anova(model)
  af$PctExp <- af$"Sum Sq"/sum(af$"Sum Sq")*100
  return(as.data.table(af, keep.rownames = "vars")[, .(vars, PctExp)])
}
