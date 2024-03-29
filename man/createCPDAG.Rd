% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{createCPDAG}
\alias{createCPDAG}
\title{Create the completed partially directed acyclic graph (CPDAG) for the input directed acyclic graph (DAG). The CPDAG represents the Markov equivalence class of the true cauasl DAG. The PC algorithms are only identifiable up to the Markov equivalence class, so assessments of causal structure recovery should be compared to the CPDAG rather than the causal DAG.}
\usage{
createCPDAG(graph)
}
\arguments{
\item{graph}{The graph object used to generate the CPDAG. Should be the ground-truth causal DAG}
}
\value{
The CPDAG corresponding to the input DAG
}
\description{
Create the completed partially directed acyclic graph (CPDAG) for the input directed acyclic graph (DAG). The CPDAG represents the Markov equivalence class of the true cauasl DAG. The PC algorithms are only identifiable up to the Markov equivalence class, so assessments of causal structure recovery should be compared to the CPDAG rather than the causal DAG.
}
\examples{
data(dag_n10000_p10)
cpdag <- rCausalMGM::createCPDAG(dag_n10000_p10)
}
