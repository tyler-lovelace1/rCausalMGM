library(bnlearn)
library(rCausalMGM)

# df0 <- read.table("data/data0.txt", header=T)
# df1 <- read.table("data/data1.txt", header=T)
# df2 <- read.table("data/data2.txt", header=T)
# df3 <- read.table("data/data3.txt", header=T)
# df4 <- read.table("data/data4.txt", header=T)

# df0 <- df0[, 1:500]

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

time_bnlearn <- function(df, sampleSize=nrow(df), filename=NULL, numTrials=5) {
    df <- df[1:sampleSize, ]
    df <- factorize_discrete(df)

    times = c()
    for (i in 1:numTrials) {
        time <- system.time( graph <- bnlearn::pc.stable(df) )[["elapsed"]]
        times <- c(times, time)
    }

    if (!is.null(filename)) {
        bnlearn::write.dot(filename, graph)
    }

    times
}

time_rCausalMGM <- function(df, sampleSize=nrow(df), filename=NULL, numTrials=5) {
    df <- df[1:sampleSize, ]
    #df <- factorize_discrete(df)

    times = c()
    for (i in 1:numTrials) {
        time <- system.time( graph <- rCausalMGM::pcStable(df) )[["elapsed"]]
        times <- c(times, time)
    }

    if (!is.null(filename)) {
        rCausalMGM::saveGraph(graph, filename)
    }

    times
}

dimensionalities <- c(25, 50, 100, 500)
sparsities <- c(2, 4, 6)
simulations <- 1:5

for (d in dimensionalities) {
    for (s in sparsities) {
        for (n in simulations) {
            dataPath <- "simulations/sim_" + d + "_" + s + "/data/data." + n + ".txt"
            dataPath
            df <- read.table(dataPath, header=T)
        }
    }
}

time_bnlearn(df0, sampleSize=1000, numTrials=1, filename="../test_results/bnlearn_test.txt")
time_rCausalMGM(df0, sampleSize=100, filename="../test_results/rCausalMGM_test.txt")
