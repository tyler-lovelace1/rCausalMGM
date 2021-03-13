#' A print override function for the graph
#'
#' @param g The graph object
#' @export
print.graph <- function(g) {
    cat("Graph generated using algorithm: ", g[["algorithm"]], "\n")
    cat("Num nodes: ", length(g[["nodes"]]), "\n")
    cat("Num edges: ", length(g[["edges"]]), "\n")
    # TODO - number of directed and instability
    cat("Num undirected edges: ", length(which(grepl("---", g[["edges"]], fixed=TRUE))), "\n")
    cat("Num directed edges: ", length(which(grepl("-->", g[["edges"]], fixed=TRUE))), "\n")
    cat("Num bidirected edges: ", length(which(grepl("---", g[["edges"]], fixed=TRUE))), "\n")
    if (!is.null(g[["stabilities"]])) {
        cat("Average instability: ", mean( 2 * g[["stabilities"]] * (1 - g[["stabilities"]]) ), "\n")
    }
}
