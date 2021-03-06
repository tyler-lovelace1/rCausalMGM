% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pcStable}
\alias{pcStable}
\title{Runs the causal algorithm PC Stable on a dataset}
\usage{
pcStable(
  df,
  maxDiscrete = 5L,
  initialGraph = NULL,
  alpha = 0.05,
  verbose = as.logical(c(FALSE))
)
}
\arguments{
\item{df}{The dataframe}

\item{maxDiscrete}{The maximum number of unique values a variable can have before being considered continuous. Defaults to 5}

\item{initialGraph}{The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.}

\item{alpha}{The p value below which results are considered significant. Defaults to 0.05.}

\item{verbose}{Whether or not to output additional information. Defaults to FALSE.}
}
\value{
The calculated search graph
}
\description{
Runs the causal algorithm PC Stable on a dataset
}
\examples{
data("data.n100.p25")
ig <- rCausalMGM::mgm(data.n100.p25)
g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
}
