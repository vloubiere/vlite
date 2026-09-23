#' Build STRING db
#' @param species either "Dm" or "Mm"
#' @param network_type The type of interactions to be used. Can be one of "full" (full functional, default) or "physical" if only physical interactions are to be considered
#' @param version databse version. default= "12.0"
#' @param input.directory The input directory where db will be saved. Default= tempdir().
#'
#' @export
stringGetDB <- function(
    species,
    network_type= "full",
    version= "12.0",
    input_directory= tempdir()
)
{
  if(input_directory!="")
    dir.create(input_directory, showWarnings = F, recursive = T)
  species.code <- switch(species, "Dm"= 7227, "Mm"= 10090, "Hs"= 9606)
  
  # Download database files if needed (see input_directory comment below)
  links <- paste0("protein.links.", network_type, ".v", version, "/", species.code, ".protein.links.v", version, ".txt.gz")
  aliases <- paste0("protein.aliases.v", version, "/", species.code, ".protein.aliases.v", version, ".txt.gz")
  info <- paste0("protein.info.v", version, "/", species.code, ".protein.info.v", version, ".txt.gz")
  for(file in c(links, aliases, info)) {
    if(!file.exists(file.path(input_directory, basename(file))))
      vlite::cmd_download(
        url = file.path("https://stringdb-downloads.org/download/", file),
        output.name = basename(file),
        output.folder = input_directory
      )
  }
  
  # Download the database
  STRINGdb::STRINGdb$new(
    version = version,
    species = switch(species, "Dm"= 7227, "Mm"= 10090, "Hs"= 9606),
    score_threshold = 0,
    network_type = network_type,
    input_directory = input_directory # Somehow changing this directory caused the download to fail for me
  )
}
