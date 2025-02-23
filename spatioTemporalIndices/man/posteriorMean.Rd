% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posteriorMean.R
\name{posteriorMean}
\alias{posteriorMean}
\title{Posterior mean. Only do the bias correction for spatial effect for catch-at-length}
\usage{
posteriorMean(run)
}
\arguments{
\item{run}{The fitted model}
}
\description{
Posterior mean. Only do the bias correction for spatial effect for catch-at-length
}
\details{
Return posterior mean of the index (using the bias-correction flag in TMB::sdreport for the catch-at-length spatial and spatio-temporal random effects only)
}
