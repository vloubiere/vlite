#' get STRING interactions Object
#'
#' Extract interactions form STRINGdb
#'
#' @param STRINGdb A STRINGdb object. See ?vl_STRING_getDB
#' @param symbols vector of Dmel gene symbols
#' @param score.cutoff If specified, only plot interactions with score>=cutoff
#' @param top.N If specified, only plot top N interactions
#' @param remove.non.connected Remove proteins that have 0 connections. Default= T
#' @param plot Should the graph be plotted?
#' @param size Vertex size. Default= 15
#' @param size2 Vertex size2 for certain shapes. Default= 15
#' @param color Vertex color. Default= "tomato"
#' @param frame.color Frame color. Default= "black"
#' @param frame.width Frame width. Default= 1
#' @param shape Shape of the vertex. Default= "circle"
#' @param label Label of the vertex. Default= symbols
#' @param label.family Default= "serif"
#' @param label.font Default= 1
#' @param label.cex Default= 1
#' @param label.dist The distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. If it is 1 then the label is displayed beside the vertex. Default= 0
#' @param label.degree It defines the position of the vertex labels, relative to the center of the vertices. It is interpreted as an angle in radian. Default= -pi/4.
#' @param label.color Default= "black"
#' @param ... Extra arguments to be passed to plot.igraph
#' @examples
#' # Build database
#' db <- vl_STRING_getDB(species= "Dm", network_type = "full", version= "11")
#'
#' # Extract interactions
#' symbols <- c("Pc", "Psc", "E(z)", "RpL10", "RpL11", "RpL12")
#' stringInteraction(STRINGdb = db,
#'                       symbols = symbols)
#' .i <- stringInteraction(STRINGdb = db,
#'                             symbols = symbols,
#'                             plot= F)
#' plot(.i)
#'
#' #Modify ploting parameters
#' print(paste0("param= ", paste0(colnames(igraph::as_data_frame(.i, what= "vertices")), collapse = ", ")))
#' igraph::V(.i)$color <- "cornflowerblue"
#' igraph::V(.i)$label.cex <- c(1,1,2,2,3,3)
#' plot(.i)
#'
#'
#' .i$V$color <- "cornflowerblue"
#' .i$V$label.cex <- seq(1, 3, length.out= 6)
#' .i$V$size <- seq(10, 60, 10)
#' plot(.i)
#'
#' @return a vl_STRING object that can easily be turned into an igraph object
#' @export
stringInteraction <- function(STRINGdb,
                              symbols,
                              score.cutoff= 400,
                              top.N= Inf,
                              remove.non.connected= T,
                              plot= T,
                              size= 15,
                              size2= 15,
                              color= "tomato",
                              frame.color= "black",
                              frame.width= 1,
                              shape= "circle",
                              label= symbols,
                              label.family= "serif",
                              label.font= 1,
                              label.cex= 1,
                              label.dist= 0,
                              label.degree= -pi/4,
                              label.color= "black",
                              ...)
{
  # Vertices
  V <- data.table(symbol= symbols,
                  size= size,
                  size2= size2,
                  color= color,
                  frame.color= frame.color,
                  frame.width= frame.width,
                  shape= shape,
                  label= label,
                  label.family= label.family,
                  label.font= label.font,
                  label.cex= label.cex,
                  label.dist= label.dist,
                  label.degree= label.degree,
                  label.color= label.color)
  V <- unique(V)

  # Get IDs (STRING put them in CAPS!)
  IDs <- as.data.table(STRINGdb$map(V, "symbol", removeUnmappedRows = TRUE))
  IDs[V[, .(symbol, caps= toupper(symbol))], symbol:= symbol, on= "symbol==caps"]

  # Full graphs
  E <- STRINGdb$get_interactions(IDs$STRING_id)
  E <- as.data.table(E)
  E[IDs, from:= i.symbol, on= "from==STRING_id"]
  E[IDs, to:= symbol, on= "to==STRING_id"]
  E <- E[, .(width= max(combined_score)), .(from, to)]

  # Checks and cutoffs
  E <- E[width>=score.cutoff]
  setorderv(E, "width", -1)
  E <- E[seq(nrow(E))<=top.N]
  E[, width:= width/999*3]

  # Keep only V and E interacting within subgroup
  if(remove.non.connected)
  {
    V <- V[symbol %in% E[, c(from, to)]]
    E <- E[from %in% V$symbol & to %in% V$symbol]
  }

  # Make Object
  obj <- list(V= V,
              E= E)
  class(obj) <- c("vl_STRING", "list")

  # PLOT
  if(plot)
    plot(obj, ...)

  # Return
  invisible(obj)
}
