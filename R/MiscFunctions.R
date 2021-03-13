#' A print override function for the graph
#'
#' @param x The graph object
#' @export
print.graph <- function(x) {
    cat("Graph generated using algorithm: ", x[["algorithm"]], "\n")
    cat("Num nodes: ", length(x[["nodes"]]), "\n")
    cat("Num edges: ", length(x[["edges"]]), "\n")
    # TODO - number of directed and instability
    cat("Num undirected edges: ", length(which(grepl("---", x[["edges"]], fixed=TRUE))), "\n")
    cat("Num directed edges: ", length(which(grepl("-->", x[["edges"]], fixed=TRUE))), "\n")
    cat("Num bidirected edges: ", length(which(grepl("<->", x[["edges"]], fixed=TRUE))), "\n")
    if (!is.null(g[["stabilities"]])) {
        cat("Average instability: ", mean( 2 * x[["stabilities"]] * (1 - x[["stabilities"]]) ), "\n")
    }
    invisible(x)
}
