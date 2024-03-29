% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmpenlkh.R
\name{pmpenlkh}
\alias{pmpenlkh}
\title{The penalized permutation likelihood approach for node-wise network selection}
\usage{
pmpenlkh(
  dat,
  group,
  lambda,
  burnin = 50000,
  nsamp = 10000,
  nintv = 200,
  niter = 50,
  eps = 0.01,
  maxcyc = 2
)
}
\arguments{
\item{dat}{a data matrix of nxnode dimension.}

\item{group}{a vector of positive integers of length ng
such that sum(group)=node.}

\item{lambda}{a vector of penalty values for network selection,
each penalty value determines a network.}

\item{burnin}{samples to be discarded in the Metropolis algorithm for sampling permutations.}

\item{nsamp}{size of the Monte Carlo sample taken for mean and variance matrix approximation.}

\item{nintv}{sampling interval in the Monte Carlo Markov chain (Metropolis).}

\item{niter}{maximum number of iterations in finding the estimator.}

\item{eps}{convergence criterion: discrepancy between successive
iterations for terminating the iteration}

\item{maxcyc}{the maximum cycle length in the candidate permutation update.}
}
\value{
Group denoting the input node clusters, a set of networks identified corresponding
to the penalty values (with within group connections coded by NA),
annotation of network connection (0=not connected etc.)
the set of networks identified with within group connections coded by 0.
}
\description{
This approach selects network structure
under the semiparametric odds ratio model
using the penalized permutation likelihood approach.
The network selection is the joint selection.
}
\details{
This method maximizes the penalized permutation likelihood
with penalty being the l^1-norm of the parameters to select
the network structure. The selection is joint selection,
NOT a node-wise selection, and is therefore more efficient statistically.
}
\examples{
\dontrun{
n=200; p=10
datmat=matrix(rnorm(n * p), ncol = p)
vargroup=c(2, 3, 1, 4)
penaltyparam=c(10, 50, 500)
pmpenlkh(dat=datmat, group=vargroup, lambda=penaltyparam)
}

}
\references{
Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.

Chen, H. Y. and Chen, J. (2020). Network selection through semiparametric odds ratio model. Manuscript.
}
