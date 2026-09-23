#' Get image metadata
#'
#' @param paths A vector of file paths for which metadata should be extracted 
#' @param tmp.dir A temp dir where to save the metadata
#' @param cleanup.cache Should the cache be cleaned up?
#'
#' @returns
#' @export
#'
#' @examples
getMetadataImage <- function(
    paths,
    tmp.dir= tempdir(),
    cleanup.cache = FALSE
) {
  # Checks ----
  stopifnot(all(grepl(".czi$|.lif$", paths)))
  stopifnot(all(sapply(paths, file.exists)))
  
  # Assemble dat ----
  dat <- data.table(path= paths)
  
  # Check if file exists ----
  dir.create(tmp.dir, showWarnings = F, recursive = T)
  tmp <- vlite::vl_cache_file(input.list = list(dat= dat), tmp.dir = tmp.dir)
  if(cleanup.cache || !file.exists(tmp)) {
    
    # For each file ----
    meta <- dat[, {
      # Import metadata
      meta <- RBioFormats::read.metadata(path, filter.metadata = T)
      stopifnot(class(meta)[1] %in% c("ImageMetadata", "ImageMetadataList"))
      if(class(meta)[1]=="ImageMetadata")
        meta <- list(meta)
      
      # Extract metadata
      meta <- lapply(
        meta,
        function(x) {
          cur <- lapply(
            x,
            function(y) {
              y <- as.data.table(y)
              if(nrow(y))
                data.table::transpose(y, keep.names = "var")
              else
                NULL
            }
          )
          rbindlist(cur)
        }
      )
      meta <- rbindlist(meta, idcol = "serie")
      setnames(meta, "V1", "value")
      
      # Compute pixel size ----
      meta[, sizeX:= as.numeric(value[var=="sizeX"]), serie]
      meta[, sizeY:= as.numeric(value[var=="sizeY"]), serie]
      meta[, sizeZ:= as.numeric(value[var=="sizeZ"]), serie]
      meta[, sizeC:= as.numeric(value[var=="sizeC"]), serie]
      if(grepl(".czi", path)) { # Apotome
        meta[, name:= gsub("(.*)-ApoTome.*", "\\1", basename(path))]
        if("ImageScaling|ImagePixelSize" %in% meta$var) { # Old software
          meta[, c("pixelX", "pixelY"):= tstrsplit(value[var=="ImageScaling|ImagePixelSize"], ",", type.convert = T)]
        } else if("Scaling|AutoScaling|CameraPixelDistance" %in% meta$var){ # Software update
          meta[, c("pixelX", "pixelY"):= tstrsplit(value[var=="Scaling|AutoScaling|CameraPixelDistance"], ",", type.convert = T)]
        }
        meta[, objective:= as.numeric(value[var=="Information|Instrument|Objective|NominalMagnification"])]
        meta[, magnification:= as.numeric(value[var=="Scaling|AutoScaling|CameraAdapterMagnification"])]
        meta[, magnification:= magnification*objective]
        meta[, pxSizeX_um:= pixelX/magnification]
        meta[, pxSizeY_um:= pixelY/magnification]
      }
      if(grepl(".lif", path)) { # Confocal
        meta[, name:= as.character(value[var=="Image name"]), serie]
        meta[, imageSizeX:= as.numeric(value[var=="Image #0|DimensionDescription #4|Length"])*1e6, serie] # Given in meters
        meta[, imageSizeY:= as.numeric(value[var=="Image #0|DimensionDescription #5|Length"])*1e6, serie] # Given in meters
        meta[, objective:= as.numeric(value[var=="Image #0|ATLConfocalSettingDefinition #0|Magnification"]), serie]
        meta[, magnification:= as.numeric(value[var=="Image #0|ATLConfocalSettingDefinition #0|Zoom"]), serie]
        meta[, pxSizeX_um:= imageSizeX/sizeX]
        meta[, pxSizeY_um:= imageSizeY/sizeY]
      }
      
      # Select important fields
      unique(meta[, .(serie, name, sizeX, sizeY, sizeZ, sizeC, objective, magnification, pxSizeX_um, pxSizeY_um)])
    }, path]
    
    # Order columns
    setcolorder(meta, "path", after = "pxSizeY_um")
    
    # Save
    saveRDS(meta, tmp)
  } else
    meta <- readRDS(tmp)
  
  print(tmp)
  # Return
  return(meta)
}
