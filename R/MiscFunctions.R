#' A print override function for the graph class
#'
#' @param x The graph object
#' @export
print.graph <- function(x, ...) {
    cat("Algorithm: ", x[["algorithm"]], "\n")
    cat("Nodes: ", length(x[["nodes"]]), "\n")
    cat("Edges: ", length(x[["edges"]]), "\n")
    if (x[["type"]] != "undirected") {
        if (x[["type"]] == "partial ancestral graph") {
            cat("  Unoriented: ", sum(grepl("o-o", x[["edges"]], fixed=TRUE)), "\n")
            cat("  Partially Oriented: ", sum(grepl("o->", x[["edges"]], fixed=TRUE)), "\n")
        } else {
            cat("  Undirected: ", sum(grepl("---", x[["edges"]], fixed=TRUE)), "\n")
        }
        cat("  Directed: ", sum(grepl("-->", x[["edges"]], fixed=TRUE)), "\n")
        cat("  Bidirected: ", sum(grepl("<->", x[["edges"]], fixed=TRUE)), "\n")
    }

    if (!is.null(x[["lambda"]])) {
        cat("lambda = {")
        for (i in 1:(length(x[["lambda"]])-1)) {
            cat(as.numeric(x[["lambda"]][i]), ", ", sep="")
        }
        cat(as.numeric(x[["lambda"]][length(x[["lambda"]])]), "}\n", sep="")
    }

    if (!is.null(x[["alpha"]])) {
        cat("alpha = ", as.numeric(x[["alpha"]]), "\n")
    }
    
    ## if (!is.null(x[["stabilities"]])) {
    ##     cat("Average instability: ", mean( 2 * x[["stabilities"]] * (1 - x[["stabilities"]]) ), "\n")
    ## }
    invisible(x)
}


#' A plot override function for the graph class
#'
#' @param x The graph object
#' @export
plot.graph <- function(x,
                       nodes = c(),
                       nodeAttr = list(shape = "ellipse"),
                       edgeAttr = list(),
                       ...) {
    ## require("Rgraphviz")

    dotParams <- list(...)
    
    nodelist <- x[["nodes"]]
    
    if (length(nodes) > 0 & !all(is.element(nodes, nodelist))) {
        stop("elements in nodes must be in the graph")
    }

    if (length(nodes)==1) {
        nodes <- c(nodes, x[["markov.blankets"]][[nodes]])
    }
    
    edges <- strsplit(sort(x[["edges"]]), " ")

    edgeL <- list()

    for (i in 1:length(nodelist)) {
        edgeL[[nodelist[i]]] <- list(edges = c(), weights = c())
    }

    edgeorient <- list(arrowhead = c(),
                       arrowtail = c())

    for (i in 1:length(edges)) {
        edgeL[[edges[[i]][1]]] <-
            list(edges = c(edgeL[[edges[[i]][1]]][["edges"]], edges[[i]][3]),
                 weights = c(edgeL[[edges[[i]][1]]][["weights"]], 1))

        edgename <- paste(edges[[i]][c(1, 3)], collapse = "~")
        edgetype <- edges[[i]][2]
        
        edgeorient[["arrowhead"]] <-
            c(edgeorient[["arrowhead"]],
              ifelse(any(edgetype == c("-->", "o->", "<->")), "open",
              ifelse(edgetype == "o-o", "odot", "none")))
        
        names(edgeorient[["arrowhead"]])[i] <- edgename

        edgeorient[["arrowtail"]] <-
            c(edgeorient[["arrowtail"]],
              ifelse(any(edgetype == c("o-o", "o->")), "odot",
              ifelse(edgetype == "<->", "open", "none")))
        
        names(edgeorient[["arrowtail"]])[i] <- edgename
    }

    rgraph <- graph::graphNEL(nodelist, edgeL, "directed")

    if (length(nodes) > 0) {
        rgraph <- graph::subGraph(nodes, rgraph)
        
        edgeorient[["arrowhead"]] <-
            edgeorient[["arrowhead"]][graph::edgeNames(rgraph)]
        
        edgeorient[["arrowtail"]] <-
            edgeorient[["arrowtail"]][graph::edgeNames(rgraph)]
    }

    graph::nodeRenderInfo(rgraph) <- nodeAttr

    edgeAttr[["arrowhead"]] <- edgeorient[["arrowhead"]]
    edgeAttr[["arrowtail"]] <- edgeorient[["arrowtail"]]
    
    graph::edgeRenderInfo(rgraph) <- edgeAttr

    if (is.null(dotParams[["main"]])) {
        dotParams[["main"]] <-
            paste0(x[["algorithm"]], 
                   ifelse(is.null(x[["lambda"]]), c(""),
                          paste0("\nlambda = {",
                                 paste(as.numeric(round(x[["lambda"]], 3)),
                                       collapse=", "), "}")),
                   ifelse(is.null(x[["alpha"]]), c(""),
                          paste0("\nalpha = ",
                                 as.numeric(x[["alpha"]]))))
    }

    graph::graph.par(list(graph = append(list(main = dotParams[["main"]]),
                                         dotParams[-which(names(dotParams) == "main")])))

    rgraph <- Rgraphviz::layoutGraph(rgraph)

    Rgraphviz::renderGraph(rgraph)
    
    ## Rgraphviz::plot(rgraph)
}

#' A table to generate a data.frame for objects from graph class. It incorporates
#' adjacency and orientation frequency if estimates of edge stability are available.
#'
#' @param graph The graph object
#' @param stabilities The stability data.frame from bootstrapping or StEPS. If NULL,
#' the stabilities entry of the graph object is used. If that is also NULL, only edge
#' interactions are returned. The default is NULL
#' @export
graphTable <- function(graph, stabilities = NULL) {

    if (is.null(stabilities)) {
        stabilities <- graph[["stabilities"]]
    }

    if (is.null(stabilities)) {
        graph.table <- data.frame(source=c(),
                                  target=c(),
                                  interaction=c())

        idx <- 1
        for (edge in graph$edges) {
            edge <- strsplit(edge, " ")[[1]]

            if (edge[2] == "---") {
                inter <- "undir"
            } else if (edge[2] == "-->") {
                inter <- "dir"
            } else if (edge[2] == "o->") {
                inter <- "ca"
            } else if (edge[2] == "<->") {
                inter <- "bidir"
            } else if (edge[2] == "o-o") {
                inter <- "cc"
            }

            tempRow <- data.frame(source=c(edge[1]),
                                  target=c(edge[3]),
                                  interaction=c(inter))
                    
            graph.table <- rbind(graph.table, tempRow)
        }
            
    } else if (graph[["type"]] == "undirected") {
        
        graph.table <- data.frame(source=c(),
                                  target=c(),
                                  interaction=c(),            
                                  adjacency.freq=c())
        
        if ((ncol(stabilities) == nrow(stabilities)) &&
            all(colnames(stabilities) == rownames(stabilities))) {

            idx <- 1
            for (edge in graph$edges) {
                edge <- strsplit(edge, " ")[[1]]
                
                tempRow <- data.frame(
                    source=c(edge[1]),
                    target=c(edge[3]),
                    interaction=c("undir"),
                    adjacency.freq=c(stabilities[edge[1], edge[3]]))
                    
                graph.table <- rbind(graph.table, tempRow)
            }
        } else {
            
            idx <- 1
            for (edge in graph$edges) {
                edge <- strsplit(edge, " ")[[1]]

                for (i in 1:nrow(stabilities)) {
                    if ((edge[1] == stabilities[i,1] &&
                         edge[3] == stabilities[i,3]) ||
                        (edge[1] == stabilities[i,3] &&
                         edge[3] == stabilities[i,1])) {

                        adj.freq <- 1-stabilities[i,"none"]
                        tempRow <- data.frame(
                            source=c(edge[1]),
                            target=c(edge[3]),
                            interaction=c("undir"),
                            adjacency.freq=c(adj.freq))
                
                        graph.table <- rbind(graph.table, tempRow)
                        
                        break
                    }
                }
            }
        }
        
    } else {
        graph.table <- data.frame(source=c(),
                                  target=c(),
                                  interaction=c(),            
                                  adjacency.freq=c(),
                                  orientation.freq=c())
        idx <- 1
        for (edge in graph$edges) {
            edge <- strsplit(edge, " ")[[1]]
            
            for (i in 1:nrow(stabilities)) {
                if (edge[1] == stabilities[i,1] &&
                    edge[3] == stabilities[i,3]) {
                    
                    adj.freq <- 1-stabilities[i,"none"]
                    if (edge[2] == "---") {
                        orient.freq <- stabilities[i,"undir"]
                        inter <- "undir"
                    } else if (edge[2] == "-->") {
                        orient.freq <- stabilities[i,"right.dir"]
                        inter <- "dir"
                    } else if (edge[2] == "o->") {
                        orient.freq <- stabilities[i,"right.partdir"]
                        inter <- "ca"
                    } else if (edge[2] == "<->") {
                        orient.freq <- stabilities[i,"bidir"]
                        inter <- "bidir"
                    } else if (edge[2] == "o-o") {
                        orient.freq <- stabilities[i,"nondir"]
                        inter <- "cc"
                    }
                    
                    tempRow <- data.frame(source=c(edge[1]),
                                          target=c(edge[3]),
                                          interaction=c(inter),
                                          adjacency.freq=c(adj.freq),
                                          orientation.freq=c(orient.freq))
                    
                    graph.table <- rbind(graph.table, tempRow)

                    break
                    
                } else if (edge[1] == stabilities[i,3] &&
                           edge[3] == stabilities[i,1]) {
                    
                    adj.freq <- 1-stabilities[i,"none"]
                    if (edge[2] == "---") {
                        orient.freq <- stabilities[i,"undir"]
                        inter <- "undir"
                    } else if (edge[2] == "-->") {
                        orient.freq <- stabilities[i,"left.dir"]
                        inter <- "dir"
                    } else if (edge[2] == "o->") {
                        orient.freq <- stabilities[i,"left.partdir"]
                        inter <- "ca"
                    } else if (edge[2] == "<->") {
                        orient.freq <- stabilities[i,"bidir"]
                        inter <- "bidir"
                    } else if (edge[2] == "o-o") {
                        orient.freq <- stabilities[i,"nondir"]
                        inter <- "cc"
                    }
                    
                    tempRow <- data.frame(source=c(edge[1]),
                                          target=c(edge[3]),
                                          interaction=c(inter),
                                          adjacency.freq=c(adj.freq),
                                          orientation.freq=c(orient.freq))
                    
                    graph.table <- rbind(graph.table, tempRow)

                    break
                }
            }
        }
    }
    
    return(graph.table)
    
}

#' A as.data.frame override function for the graph class
#'
#' @param x The graph object
#' @export
as.data.frame.graph <- function(x, ...) {
    return(graphTable(x))
}
