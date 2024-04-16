#' A print override function for the graph class
#'
#' @param x The graph object
#' @param ... Additional print arguments
#' @export
print.graph <- function(x, ...) {
    cat("Algorithm: ", x[["algorithm"]], "\n")
    cat("Nodes: ", length(x[["nodes"]]), "\n")
    cat("Edges: ", length(x[["edges"]]), "\n")
    if (x[["type"]] != "undirected") {
        if (x[["type"]] == "partial ancestral graph") {
            cat("  Unoriented: ", sum(grepl("o-o", x[["edges"]], fixed=TRUE)), "\n")
            cat("  Partially Oriented: ", sum(grepl("o->", x[["edges"]], fixed=TRUE)), "\n")
            cat("  Directed: ", sum(grepl("-->", x[["edges"]], fixed=TRUE)), "\n")
            cat("  Bidirected: ", sum(grepl("<->", x[["edges"]], fixed=TRUE)), "\n")
        } else {
            cat("  Directed: ", sum(grepl("-->", x[["edges"]], fixed=TRUE)), "\n")
            cat("  Undirected: ", sum(grepl("---", x[["edges"]], fixed=TRUE)), "\n")
        }
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

    if (!is.null(x[["penalty"]])) {
        cat("penalty = ", as.numeric(x[["penalty"]]), "\n")
    }
    
    if (!is.null(attr(x, "Score"))) {
        cat("score = ", as.numeric(attr(x, "Score")), "\n")
    }
    
    invisible(x)
}

#' A print override function for the graphCV class
#'
#' @param x The graphCV object
#' @param ... Additional print arguments
#' @export
print.graphCV <- function(x, ...) {
    mgmFlag <- !is.null(x$lambdas)
    cat("Minimum cross-validation log(pseudo-likelihood) graph:\n")
    cat("  Index: ")
    cat(x$idx.min)
    cat("\n  Average Markov Blanket Size: ")
    cat(mean(sapply(x$graph.min$markov.blankets, length)))
    cat("\n\n")
    print(x$graph.min)
    cat("\n")
    cat("One standard error above the minimum cross-validation log(pseudo-likelihood) graph:\n")
    cat("  Index: ")
    cat(x$idx.1se)
    cat("\n  Average Markov Blanket Size: ")
    cat(mean(sapply(x$graph.1se$markov.blankets, length)))
    cat("\n\n")
    print(x$graph.1se)
    invisible(x)
}


#' A print override function for the graphPath class
#'
#' @param x The graphPath object
#' @param ... Additional print arguments
#' @export
print.graphPath <- function(x, ...) {
    mgmFlag <- !is.null(x$lambdas)
    cat("Minimum BIC graph:\n")
    cat("  Index: ")
    cat(which.min(x$BIC))
    cat("\n\n")
    print(x$graph.bic)
    cat("\n")
    cat("Minimum AIC graph:\n")
    cat("  Index: ")
    cat(which.min(x$AIC))
    cat("\n\n")
    print(x$graph.aic)
    invisible(x)
}


#' A print override function for the graphSTEPS class
#'
#' @param x The graphSTEPS object
#' @param ... Additional print arguments
#' @export
print.graphSTEPS <- function(x, ...) {
    mgmFlag <- !is.null(x$lambdas)
    
    cat("StEPS selected graph:\n\n")
    print(x$graph.steps)
    cat("\n\n")

    cat("StARS selected graph:\n\n")
    print(x$graph.stars)

    invisible(x)
}


#' A print override function for the graphSTARS class
#'
#' @param x The graphSTARS object
#' @param ... Additional print arguments
#' @export
print.graphSTARS <- function(x, ...) {
    cat("StARS selected graph:\n\n")
    print(x$graph)
    invisible(x)
}

#' A print override function for the knowledge class
#'
#' @param x The knowledge object
#' @param ... Additional print arguments
#' @export
print.knowledge <- function(x, ...) {
    cat("Prior Knowledge:\n")
    if (length(x$tiers) > 0) {
        cat("  Tiers:\n")
        for (t in 1:length(x$tiers)) {
            cat(paste0("    ", t,
                       ifelse(x$forbiddenWithinTier[t], '*: ', ': '),
                       paste(x$tiers[[t]], collapse=" "), "\n"))
        }
    }
    if (length(x$forbidden) > 0 ) {
        cat("  Forbidden:\n")
        for (pair in x$forbidden) {
            cat(paste0("    ", pair[1], " --> ", pair[2], "\n"))
        }
    }
    if (length(x$required) > 0 ) {
        cat("  Required:\n")
        for (pair in x$required) {
            cat(paste0("    ", pair[1], " --> ", pair[2], "\n"))
        }
    }
    invisible(x)
}

#' A plot override function for the graph class
#'
#' @param x The graph object
#' @param nodes A subset of nodes in the graph to plot. If only a
#'     single node is supplied, then that node and its Markov blanket
#'     will be plotted.
#' @param nodeAttr A list of options to modify graph nodes
#'     (e.g. fontsize).
#' @param edgeAttr A list of options to modify graph edges.
#' @param ... Additional plot arguments
#' @export
plot.graph <- function(x,
                       nodes = c(),
                       nodeAttr = list(),
                       edgeAttr = list(),
                       ...) {
    ## if (!require(Rgraphviz, quietly = TRUE))
    ##     install.packages("Rgraphviz")

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

    if (length(nodeAttr) > 0) {
        graph::nodeRenderInfo(rgraph) <- nodeAttr
    }
    
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
                                 as.numeric(round(x[["alpha"]], 4)))))
    }

    graph::graph.par(list(graph = append(list(main = dotParams[["main"]]),
                                         dotParams[-which(names(dotParams) == "main")])))

    rgraph <- Rgraphviz::layoutGraph(rgraph)

    Rgraphviz::renderGraph(rgraph)
    
}

#' A plot override function for the graphCV class
#'
#' @param x The graph object
#' @param ... Additional plot arguments
#' @export
plot.graphCV <- function(x, ...) {
    mgmFlag <- x$graph.min$type=="undirected"
    
    if (mgmFlag) {
        log10params <- log10(x$lambdas)
        minParam <- log10(x$lambdas[x$idx.min])
        seParam <- log10(x$lambdas[x$idx.1se])
    } else {
        log10params <- x$size
        minParam <- x$size[x$idx.min]
        seParam <- x$size[x$idx.1se]
    }
    
    llMeans <- x$mean
    llSe <- x$se

    upperVal <- stats::quantile(llMeans+llSe, probs=0.75)
    lowerVal <- stats::quantile(llMeans-llSe, probs=0.25)

    iqr <- upperVal - lowerVal
    ll.lims <- c(lowerVal - 1.5 * iqr,
                 upperVal + 1.5 * iqr)

    plot(x=log10params, y=llMeans, col='red', pch=19,
         xlab=ifelse(mgmFlag, expression(log10(lambda)), "Average Markov Blanket Size"),
         ylab="-log(Pseudo-Likelihood)", ylim=ll.lims)
    
    graphics::arrows(x0=log10params, x1=log10params, y0=llMeans-llSe, code=3, angle=90,
                     length=0.05, y1=llMeans+llSe, col='darkgray')
    
    graphics::abline(v=minParam,
                     col='black', lty=3, lw=2)
    
    graphics::abline(v=seParam,
                     col='black', lty=3, lw=2)
}

#' A plot override function for the graphCV class
#'
#' @param x The graph object
#' @param ... Additional plot arguments
#' @export
plot.graphPath <- function(x, ...) {
    mgmFlag <- !is.null(x$lambdas)

    if (mgmFlag) {
        log10params <- log10(x$lambdas)
    } else {
        log10params <- log10(x$alphas)
    }
    
    score.range <- (max(x$BIC) - min(x$AIC)) / (2 * x$n)
    score.lims <- c(min(x$AIC) / (2 * x$n) - 0.025 * score.range,
                    max(x$BIC) / (2 * x$n) + 0.025 * score.range)

    plot(x=log10params, y=x$BIC / (2 * x$n), col='red', pch=19,
         xlab=ifelse(mgmFlag, expression(log10(lambda)), expression(log10(alpha))),
         ylab="Sample Averaged Score", ylim=score.lims)

    graphics::points(x=log10params, y=x$AIC / (2 * x$n), col='blue', pch=19)

    graphics::legend(x = "bottomright", title="Scores", 
                     legend=c("AIC", "BIC"), 
                     col = c("blue","red"),
                     pch=19, cex=0.7)

    graphics::abline(v=ifelse(mgmFlag,
                              log10(x$lambdas[which.min(x$AIC)]),
                              log10(x$alphas[which.min(x$AIC)])),
                     col='blue', lty=3, lw=2)
    
    graphics::abline(v=ifelse(mgmFlag,
                              log10(x$lambdas[which.min(x$BIC)]),
                              log10(x$alphas[which.min(x$BIC)])),
                     col='red', lty=3, lw=2)
}

#' A plot override function for the graphSTEPS class
#'
#' @param x The graph object
#' @param ... Additional plot arguments
#' @export
plot.graphSTEPS <- function(x, ...) {
    mgmFlag <- !is.null(x$lambdas)

    gamma <- x$gamma

    if (mgmFlag) {
        log10params <- log10(x$lambdas)
    } else {
        log10params <- log10(x$alphas)
    }

    plot(x=log10params, y=x$instability[,4], col='black', pch=19,
         xlab=ifelse(mgmFlag, expression(log10(lambda)), expression(log10(alpha))),
         ylab="Edge instability across subsamples", ylim=c(0,0.5))

    graphics::points(x=log10params, y=x$instability[,1], col='red', pch=18)
    graphics::points(x=log10params, y=x$instability[,2], col='blue', pch=17)
    graphics::points(x=log10params, y=x$instability[,3], col='purple', pch=15)

    graphics::abline(h=gamma, lty=2, col='gray', lw=2)
    
    graphics::abline(v=log10(x$graph.steps$lambda[1]),  col='red',    lty=2, lw=2)
    graphics::abline(v=log10(x$graph.steps$lambda[2]),  col='blue',   lty=2, lw=2)
    graphics::abline(v=log10(x$graph.steps$lambda[3]),  col='purple', lty=2, lw=2)
    graphics::abline(v=log10(x$graph.stars$lambda[1]),  col='black',  lty=3, lw=2)

    graphics::legend(x = "topleft", title="Edge Type", 
                     legend = c("All", "C-C", "C-D", "D-D"), 
                     col = c("black","red", "blue", "purple"),
                     pch = c(19, 18, 17, 15), cex=0.7)
}

#' A plot override function for the graphSTARS class
#'
#' @param x The graph object
#' @param ... Additional plot arguments
#' @export
plot.graphSTARS <- function(x, ...) {
    gamma <- x$gamma
    log10params <- log10(x$alphas)
    
    plot(x=log10params, y=x$instability[,1], col='black', pch=19,
         xlab=expression(log10(alpha)),
         ylab="Edge instability across subsamples", ylim=c(0,0.1))

    graphics::abline(h=gamma, lty=2, col='gray', lw=2)

    graphics::abline(v=log10(x$graph$alpha),  col='black',  lty=3, lw=2)
}


#' A function to generate a data.frame for objects from graph
#' class. It incorporates adjacency and orientation frequency if
#' estimates of edge stability are available.
#'
#' @param graph The graph object
#' 
#' @param stabilities The stability data.frame from bootstrapping or
#'     StEPS. If NULL, the stabilities entry of the graph object is
#'     used. If that is also NULL, only edge interactions are
#'     returned. The default is NULL
#'
#' @return A data.frame containing source, target, and interaction
#'     columns for each edge in the graph. If stabilities are
#'     available, then the adjFrequency and orientation
#'     frequencies (if applicable) are returned for each edge.
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
                                  adjFreq=c())
        
        if ((ncol(stabilities) == nrow(stabilities)) &&
            all(colnames(stabilities) == rownames(stabilities))) {

            idx <- 1
            for (edge in graph$edges) {
                edge <- strsplit(edge, " ")[[1]]
                
                tempRow <- data.frame(
                    source=c(edge[1]),
                    target=c(edge[3]),
                    interaction=c("undir"),
                    adjFreq=c(stabilities[edge[1], edge[3]]))
                
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
                            adjFreq=c(adj.freq))
                        
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
                                  adjFreq=c(),
                                  orientFreq=c())
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
                                          adjFreq=c(adj.freq),
                                          orientFreq=c(orient.freq))
                    
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
                                          adjFreq=c(adj.freq),
                                          orientFreq=c(orient.freq))
                    
                    graph.table <- rbind(graph.table, tempRow)

                    break
                }
            }
        }
    }
    
    return(graph.table)
    
}


#' A function to create a prior knowledge object for use with causal
#' discovery algorithms
#'
#' @param tiers A list containing ordered vectors of variables where
#'     variables in tier t can only be ancestors of variables in tiers
#'     t+1 ... T and descendants of variables in tiers (1 .. t-1). If
#'     tiers are used, all variables must be in a tier, and no
#'     variable can be in multiple tiers.
#' @param forbiddenWithinTier A vector of logical values indicating
#'     whether edges are allowed between variables in a given
#'     tier. The value is NULL by default, which results in
#'     forbiddenWithinTier being set to FALSE for each tier.
#' @param forbidden A list containing vectors of node pairs that
#'     forbid a specific directed edge. For example, to forbid 
#'     A --> B, add c("A", "B") to forbidden.
#' @param required A list containing vectors of node pairs that
#'     require the presence of a specific directed edge. For example,
#'     to require B --> A, add c("B", "A") to required.
#' @return A knowledge object that can be passed to causal discovery
#'     algorithms.
#' @export
createKnowledge <- function(tiers = list(), forbiddenWithinTier=NULL,
                            forbidden=list(), required=list()) {
    if (is.null(forbiddenWithinTier)) {
        if (length(tiers) != 0) {
            forbiddenWithinTier <- rep(FALSE, length(tiers))
        } else {
            forbiddenWithinTier <- as.logical(c())
        }
    }
    knowledge <- list(tiers=tiers, forbiddenWithinTier=forbiddenWithinTier,
                      forbidden=forbidden, required=required)
    class(knowledge) <- 'knowledge'
    return(knowledge)
}


#' A function to simulate a random forward DAG from a SEM model.
#'
#' @param n The sample size of the generated dataset. The default is
#'     1000.
#' @param p The number of features in the generated dataset. The
#'     default is 50.
#' @param discFrac The fraction of variables in the dataset that are
#'     discrete. The default is 0.5.
#' @param deg The average graph degree for the simulated graph. The
#'     default is 3.
#' @param coefMin The lower bound on the magnitude of the effect
#'     size. The default is 0.5.
#' @param coefMax The upper bound on the magnitude of the effect
#'     size. The default is 1.5.
#' @param noiseMin The lower bound on the standard deviation of the
#'     Gaussian noise for continuous variables. The default is 1.
#' @param noiseMax The upper bound on the standard deviation of the
#'     Gaussian noise for continuous variables. The default is 2.
#' @param seed The random seed for generating the simulated DAG. The
#'     default is NULL.
#' @return A list containing the simulated dataset and the
#'     corresponding ground truth causal DAG.
#' @examples
#' sim <- simRandomDAG(200, 25)
#' print(sim$graph)
#' print(sim$data[1:6,])
#' @export
simRandomDAG <- function(n=1000, p=50, discFrac=0.5, deg=3,
                         coefMin=0.5, coefMax=1.5, noiseMin=1, noiseMax=2,
                         seed=NULL) {

    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    numCat <- floor(p * discFrac)
    numCont <- p - numCat
    nodes <- c()
    if (numCont > 0) {
        nodes <- c(nodes, paste0('X',1:numCont))
    }
    if (numCat>0) {
        nodes <- c(nodes, paste0('Y',1:numCat))
    }

    permNodes <- sample(nodes)

    numEdges <- floor(deg * p / 2)

    edgeIdx <- sample(1:(p*(p-1)/2), numEdges)

    adjMat <- matrix(0, p, p)

    adjMat[lower.tri(adjMat)][edgeIdx] <- 1

    logsumexp <- function(x) {
        b <- max(x)
        x <- x-b
        return(log(sum(exp(x)))+b)
    }

    softmax <- function(x) {
        return(exp(x-logsumexp(x)))
    }

    data <- data.frame()

    idx <- 1
    for (node in permNodes) {
        pa <- permNodes[adjMat[idx,]==1]
        if (length(pa) > 0) {
            f <- stats::as.formula(paste('~', paste(pa, collapse=' + '), sep=' + '))
            mod.mat <- as.matrix(stats::model.matrix(f, data)[,-1])
            if (ncol(mod.mat)==1) {
                colnames(mod.mat) <- pa
            }
        }
        
        if (grepl('X', node)) {
            if (length(pa)==0) {
                val <- stats::rnorm(n)
            } else {
                betaScale <- stats::runif(length(pa), coefMin, coefMax)
                names(betaScale) <- pa
                beta <- matrix(stats::runif(ncol(mod.mat), -1, 1), 1)
                colnames(beta) <- colnames(mod.mat)
                for (paNode in pa) {
                    paIdx <- grep(paNode, colnames(beta))
                    if (grepl('X', paNode)) {
                        beta[,paIdx] <- sign(beta[,paIdx]) * betaScale[paNode]
                    } else {
                        beta[,paIdx] <- beta[,paIdx] - mean(beta[,paIdx])
                        beta[,paIdx] <- betaScale[paNode] * beta[,paIdx] / sqrt(sum(beta[,paIdx]^2))
                    }
                }
                
                pred <- mod.mat %*% t(beta)

                val <- as.vector(scale(pred + stats::rnorm(n, sd=stats::runif(1, noiseMin, noiseMax))))
            }
        } else if (grepl('Y', node)) {
            if (length(pa)==0) {
                logprobs <- matrix(rep(0, 3*n), n)
                val <- factor(
                    apply(logprobs, 1, function(x) sample(c('A','B','C'), 1, prob=softmax(x))),
                    levels=c('A','B','C')
                )
            } else {
                betaScale <- stats::runif(length(pa), coefMin, coefMax)
                names(betaScale) <- pa
                
                beta <- matrix(stats::runif(3 * ncol(mod.mat), -1, 1), 3)
                colnames(beta) <- colnames(mod.mat)
                for (paNode in pa) {
                    paIdx <- grep(paNode, colnames(beta))
                    beta[,paIdx] <- beta[,paIdx] - mean(beta[,paIdx])
                    beta[,paIdx] <- betaScale[paNode] * beta[,paIdx] / sqrt(sum(beta[,paIdx]^2))
                }

                logprobs <- mod.mat %*% t(beta)

                val <- factor(
                    apply(logprobs, 1, function(x) sample(c('A','B','C'), 1, prob=softmax(x))),
                    levels=c('A','B','C')
                )
            }
        }
        if (idx==1) {
            data <- data.frame(val)
            colnames(data) <- node
        } else {
            data[,node] <- val
        }
        idx <- idx + 1
    }
    
    graph <- adjMat2Graph(adjMat, permNodes, directed=T)

    graph$algorithm <- "Ground Truth"
    graph$type <- "directed acyclic graph"

    result <- list(data=data[,nodes], graph=graph)
    return(result)
}
