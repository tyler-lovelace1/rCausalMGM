% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{prMetricsOrientation}
\alias{prMetricsOrientation}
\title{Orientation Precision-Recall Metrics}
\usage{
prMetricsOrientation(estimate, groundTruth, groundTruthDAG = NULL)
}
\arguments{
\item{estimate}{An estimated graph object}

\item{groundTruth}{A ground truth graph object of the same type as the estimated graph object}

\item{groundTruthDAG}{A ground truth graph object containing the true causal DAG. Only necessary for calculating the or precision, recall, F1, and MCC for partial ancestral graphs (PAGs)}
}
\value{
The orientation precision, recall, F1, and MCC, between the two graph objects
}
\description{
Calculate the orientation precision, recall, F1, and Matthew's Correlation Coefficient (MCC) between an estimated and ground truth graph.
}
\examples{
data("train_n10000_p10")
sim <- simRandomDAG(200, 25)
g <- pcStable(sim$data)
prMetricsOrientation(g, cpdag(sim$graph))
}
