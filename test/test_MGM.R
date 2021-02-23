library(rCausalMGM)

# df <- read.table("data/data5.csv", header=T) # Tiny dataset
df <- read.table("data/data.n100.p25.txt", header=T)
# df <- df[-which(df[,"X128"]==3),]
# df <- df[-which(df[,"X274"]==1),]

lam.max <- max(abs(cor(df) - diag(ncol(df))))
lam.min <- 0.2*lam.max

lams <- exp(seq(log(lam.max), log(lam.min), -(log(lam.max) - log(lam.min))/14))
# lams
# rCausalMGM::runTests(df)

# g <- rCausalMGM::loadGraph("../cpc_out_censored.txt")
# rCausalMGM::saveGraph(g, "../r_graph_out.txt")
# rCausalMGM::pcStable(df, verbose=TRUE)

g <- rCausalMGM::steps(df, lambda = lams, computeStabs=TRUE, verbose=TRUE)
print(g)
rCausalMGM::saveGraph(g, "../testStabGraph.txt")
# rCausalMGM::saveGraph(g, '../mgm_local.txt')
# g.cpc <- rCausalMGM::cpcStable(df, initialGraph = g, verbose = TRUE)
# print(g.cpc)
# rCausalMGM::saveGraph(g.cpc, '../cpc_local.txt')
# g.pcm <- rCausalMGM::pcMax(df, initialGraph = g, verbose = TRUE)
# print(g.pcm)
# rCausalMGM::saveGraph(g.pcm, '../pcm_local.txt')

# rCausalMGM::pcMax(df, lambda = c(0.05, 0.05, 0.05), verbose = TRUE)

