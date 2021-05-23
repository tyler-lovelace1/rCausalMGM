library(rCausalMGM)
# library(foreach)
data(data.n1000.p100)
# bs.graphs <- foreach(i = 1:100) %do% {
#     idxs <- sample(1:nrow(data.n100.p25), nrow(data.n100.p25), replace=T)
#     ig <- mgm(data.n100.p25[idxs,], lambda = c(0.15, 0.15, 0.15))
#     g <- pcMax(data.n100.p25[idxs,], initialGraph = ig)
# }

ig <- mgm(data.n1000.p100)
g <- pcMax(data.n1000.p100, initialGraph=ig, verbose = TRUE)

g <- pcMax(data.n1000.p100, initialGraph=ig, threads = 2, verbose = TRUE)


# saveGraph(g, "../fci.txt")

# fg <- loadGraph("../fci.txt")
# print(fg)

# fg[["markov.blankets"]]

# lam.max <- max(abs(cor(df) - diag(ncol(df))))
# lam.min <- 0.2*lam.max

# lams <- exp(seq(log(lam.max), log(lam.min), -(log(lam.max) - log(lam.min))/14))

# g <- rCausalMGM::steps(df, lambda = lams, leaveOneOut=TRUE, computeStabs=TRUE, verbose=TRUE)
# print(g)
# rCausalMGM::saveGraph(g, "../testStabGraph.txt")

# rCausalMGM::pcStable(df, initialGraph=g)
# rCausalMGM::cpcStable(df, initialGraph=g)
# rCausalMGM::pcMax(df, initialGraph=g)

# g <- rCausalMGM::loadGraph("graph/graph.n1000.p100.txt")
# g[["markov.blankets"]]


