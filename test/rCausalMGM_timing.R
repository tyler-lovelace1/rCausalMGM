library(bnlearn)
library(rCausalMGM)
library(parallel)

# df0 <- read.table("data/data0.txt", header=T)
# df1 <- read.table("data/data1.txt", header=T)
# df2 <- read.table("data/data2.txt", header=T)
# df3 <- read.table("data/data3.txt", header=T)
# df4 <- read.table("data/data4.txt", header=T)

# df0 <- df0[, 1:500]

options(warn=-1)

NCPUs = 32

# In bn_learn, integer columns need to be converted to factors
factorize_discrete <- function(df, maxDiscrete = 5) {
    discrete_vars = c()
    for(var in colnames(df)) {
        if (length(unique(df[[var]])) <= maxDiscrete) {
            discrete_vars <- c(discrete_vars, var)
            # print(var)
            # print(length(unique(df_small[[var]])))
        }
    }
    df[,discrete_vars] <- lapply(df[,discrete_vars], as.factor)

    df
}


timing <- function(n, p, degree, cluster=NULL, vs = 1:5) {

    


    for (v in vs) {
        dataPath <- "data/data.n100.p25.txt"

        # path.bnlearn   <- paste("simulations/sim_", p, "_", degree, "/graph/bnlearn.",   n, ".", v, ".txt", sep="")
        # path.pcStable   <- paste("simulations/sim_", p, "_", degree, "/graph/pcStable.",   n, ".", v, ".txt", sep="")
        # path.cpcStable  <- paste("simulations/sim_", p, "_", degree, "/graph/cpcStable.",  n, ".", v, ".txt", sep="")
        # path.pcMax      <- paste("simulations/sim_", p, "_", degree, "/graph/pcMax.",      n, ".", v, ".txt", sep="")
        # path.pc50       <- paste("simulations/sim_", p, "_", degree, "/graph/pc50.",       n, ".", v, ".txt", sep="")
        # path.fciStable  <- paste("simulations/sim_", p, "_", degree, "/graph/fciStable.",  n, ".", v, ".txt", sep="")
        # path.cfciStable <- paste("simulations/sim_", p, "_", degree, "/graph/cfciStable.", n, ".", v, ".txt", sep="")
        # path.fciMax     <- paste("simulations/sim_", p, "_", degree, "/graph/fciMax.",     n, ".", v, ".txt", sep="")

        df <- read.table(dataPath, header=T)

        time <- round(system.time( { ig <- rCausalMGM::mgm(df); graph <- rCausalMGM::pcStable(df, initialGraph=ig) } )[["elapsed"]])
        # rCausalMGM::saveGraph(graph, path.pcStable)
        cat(paste("pcStable", n, p, degree, NCPUs, v, time, sep=","))
        cat("\n")

        time <- round(system.time( { ig <- rCausalMGM::mgm(df); graph <- rCausalMGM::cpcStable(df, initialGraph=ig) } )[["elapsed"]])
        # rCausalMGM::saveGraph(graph, path.cpcStable)
        cat(paste("cpcStable", n, p, degree, NCPUs, v, time, sep=","))
        cat("\n")

        time <- round(system.time( { ig <- rCausalMGM::mgm(df); graph <- rCausalMGM::pcMax(df, initialGraph=ig) } )[["elapsed"]])
        # rCausalMGM::saveGraph(graph, path.pcMax)
        cat(paste("pcMax", n, p, degree, NCPUs, v, time, sep=","))
        cat("\n")

        time <- round(system.time( { ig <- rCausalMGM::mgm(df); graph <- rCausalMGM::pc50(df, initialGraph=ig) } )[["elapsed"]])
        # rCausalMGM::saveGraph(graph, path.pc50)
        cat(paste("pc50", n, p, degree, NCPUs, v, time, sep=","))
        cat("\n")

        time <- round(system.time( { ig <- rCausalMGM::mgm(df); graph <- rCausalMGM::fciStable(df, initialGraph=ig) } )[["elapsed"]])
        # rCausalMGM::saveGraph(graph, path.fciStable)
        cat(paste("fciStable", n, p, degree, NCPUs, v, time, sep=","))
        cat("\n")

        time <- round(system.time( { ig <- rCausalMGM::mgm(df); graph <- rCausalMGM::cfci(df, initialGraph=ig) } )[["elapsed"]])
        # rCausalMGM::saveGraph(graph, path.cfciStable)
        cat(paste("cfciStable", n, p, degree, NCPUs, v, time, sep=","))
        cat("\n")

        time <- round(system.time( { ig <- rCausalMGM::mgm(df); graph <- rCausalMGM::fciMax(df, initialGraph=ig) } )[["elapsed"]])
        # rCausalMGM::saveGraph(graph, path.fciMax)
        cat(paste("fciMax", n, p, degree, NCPUs, v, time, sep=","))
        cat("\n")

        # df <- factorize_discrete(df)
        # time <- round(system.time( graph <- bnlearn::pc.stable(dfuster) )[["elapsed"]])
        # time.bnlearn <- c(time.bnlearn, time)
        # bnlearn::write.dot(path.bnlearn, graph)

        # gc()

    }
    
}

# cl <- makeCluster(32)

NCPUs = 32

# cat(c("Algorithm", "n", "p", "degree", "nCPUs", "v", "runtime (s)"), sep=",")
# cat("\n")

timing(100, 100, 4)
timing(250, 100, 4)
timing(500, 100, 4)
timing(1000, 100, 4)
timing(5000, 100, 4)
timing(10000, 100, 4)

timing(1000, 25, 4)
timing(1000, 50, 4)
timing(1000, 500, 4)

timing(1000, 100, 2)
timing(1000, 100, 6)

# stopCluster(cl)
