% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cdimpute.R
\name{cdimpute}
\alias{cdimpute}
\title{Impute missing values by random draw from the conditional distributions}
\usage{
cdimpute(dat = dat, miscode = c(-9), method = "sp", niter = 10, nimpute = 5)
}
\arguments{
\item{dat}{a data matrix with possible missing values}

\item{miscode}{a set of values used for denoting missing values}

\item{method}{specific method of estimation, can be
method='pw' for pairwise likelihood approach
method='sp' for semiparametric likelihood approach
method='pm' for permutation likelihood approach}

\item{niter}{the number of iterations to produce the final imputed data set}

\item{nimpute}{the number of imputed copies for each missing value}
}
\value{
\enumerate{
\item One imputed data set
}
}
\description{
This is a repeated imputation by the conditional distribution
}
\details{
This function imputes missing data by random draw from
the conditional distribution of the variable conditional on other variables.
in each iteration, the conditional model is estimated by
the semiparametric OR model using the semiparametric likelihood approach
based on the current imputed data
}
\examples{
\dontrun{
 cdimpute(dat=bone,miss=c(-9))
}

}
\references{
Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
}
