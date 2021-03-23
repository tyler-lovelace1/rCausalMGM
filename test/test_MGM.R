library(rCausalMGM)

df <- read.table("data/data.n100.p25.txt", header=T)
g <- rCausalMGM::mgm(df)
print(g)

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


