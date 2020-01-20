library(rCausalMGM)

df <- data.frame(
  C1 = c(-1.2, 0.3, 37, -2.2, 6.6, 7.96),
  C2 = c(1.2, -0.3, 17, 6.2, 1, 0),
  D1 = c(0, 1, 0, 2, 2, 0),
  D2 = c("yes", "no", "no", "yes", "yes", "no")
)

df

cmgmds1 <- rCausalMGMData(df)

df <- data.frame(
  C1 = c(-1.2, 0.3, 3.7, -2.2, 6.6, 7.96),
  C2 = c(1.2, -0.3, -2.0, 6.2, 1, -7.3),
  D1 = c(0.0, 1.0, 0.0, 2.0, 1.0, 2.0),
  D2 = c("yes", "no", "no", "yes", "yes", "maybe")
)

df

cmgmds2 <- rCausalMGMData(df, 3)

cmgmds1

cmgmds2
