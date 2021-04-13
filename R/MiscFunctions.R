#' A print override function for the graph
#'
#' @param x The graph object
#' @export
print.graph <- function(x) {
    cat("Algorithm: ", x[["algorithm"]], "\n")
    cat("Nodes: ", length(x[["nodes"]]), "\n")
    cat("Edges: ", length(x[["edges"]]), "\n")
    if (x[["type"]] == "partial ancestral graph") {
        cat("  Nondirected: ", length(which(grepl("o-o", x[["edges"]], fixed=TRUE))), "\n")
        cat("  Partially Oriented: ", length(which(grepl("o->", x[["edges"]], fixed=TRUE))), "\n")
    } else {
        cat("  Undirected: ", length(which(grepl("---", x[["edges"]], fixed=TRUE))), "\n")
    }
    cat("  Directed: ", length(which(grepl("-->", x[["edges"]], fixed=TRUE))), "\n")
    cat("  Bidirected: ", length(which(grepl("<->", x[["edges"]], fixed=TRUE))), "\n")
    
    if (!is.null(x[["stabilities"]])) {
        cat("Average instability: ", mean( 2 * x[["stabilities"]] * (1 - x[["stabilities"]]) ), "\n")
    }
    invisible(x)
}
