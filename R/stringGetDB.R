#' Build STRING db
#' @param species either "Dm" or "Mm"
#' @param network_type The type of interactions to be used. Can be one of "full" (full functional, default) or "physical" if only physical interactions are to be considered
#' @param version databse version. default= "10"
#'
#' @export
stringGetDB <- function(species,
                        network_type= "full",
                        version= "11.0")
{
  if(!network_type %in% c("full", "physical"))
    stop("Network_type has to be one of 'full' (full functional annot) or 'physical' (physical interactions)")

  STRINGdb::STRINGdb$new(version = version,
                         species = switch(species, "Dm"= 7227, "Mm"= 10090, "Hs"= 9606),
                         score_threshold = 0,
                         network_type = network_type,
                         input_directory = "")
}
