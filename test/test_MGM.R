library(rCausalMGM)

df <- read.table("rCausalMGM/data/data.n1000.p100.txt", header=T) 

#df <- read.table("rCausalMGM/data/data.n100.p25.txt", header=T)
g <- rCausalMGM::mgm(df)
# df <- read.table("rCausalMGM/data/data0.txt", header=T)
# df <- df[-which(df[,"X128"]==3),]
# df <- df[-which(df[,"X274"]==1),]

lam.max <- max(abs(cor(df) - diag(ncol(df))))
lam.min <- 0.1*lam.max

lams <- exp(seq(log(lam.min), log(lam.max), (log(lam.max) - log(lam.min))/14))

# rCausalMGM::pcStable(df, verbose=TRUE)

# print(df)
# g <- rCausalMGM::steps(df, lambda = lams, verbose=TRUE)
g <- rCausalMGM::mgm(df, verbose=TRUE)
print(g)
rCausalMGM::saveGraph(g, '../mgm_local.txt')
#g.cpc <- rCausalMGM::cpcStable(df, initialGraph = g, verbose = TRUE)
#print(g.cpc)
#rCausalMGM::saveGraph(g.cpc, '../cpc_local.txt')
#g.pcm <- rCausalMGM::pcMax(df, initialGraph = g, verbose = TRUE)
#print(g.pcm)
#rCausalMGM::saveGraph(g.pcm, '../pcm_local.txt')

g.fci <- rCausalMGM::fciStable(df, initialGraph = g, verbose = TRUE)
rCausalMGM::saveGraph(g.fci, '../fci_local.txt')

# rCausalMGM::pcMax(df, lambda = c(0.05, 0.05, 0.05), verbose = TRUE)

