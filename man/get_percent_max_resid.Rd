% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PercentMaximumDifference.R
\name{get_percent_max_resid}
\alias{get_percent_max_resid}
\title{get_percent_max_resid}
\usage{
get_percent_max_resid(observed, expected)
}
\arguments{
\item{observed}{observed matrix (i.e. contingency table)}

\item{expected}{expected matrix under the null hypothesis that there are no differences}
}
\value{
the percentage of asymmetry relative to the maximum possible asymmetry.
}
\description{
\code{get_percent_max_resid} Takes in the observed and expected matrix (as with \code{chisq.test}) and uses it to calculate the PMD.
}
\examples{
   dummy_mat <- matrix(c(500, 250, 100, 1000, 500, 200),nrow = 3, ncol=2)
   get_percent_max_resid(dummy_mat, get_expected(dummy_mat))
}
