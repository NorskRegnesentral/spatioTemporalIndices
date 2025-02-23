% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{includeIntPoints}
\alias{includeIntPoints}
\title{includeIntPoints}
\usage{
includeIntPoints(data, conf_l, confPred, gamSetup_depth)
}
\arguments{
\item{data}{data in model}

\item{conf_l}{configurations in model}

\item{confPred}{prediction configurations}

\item{gamSetup_depth}{gam setup for depth, needed to set up the spline for integration points}
}
\value{
This function includes integration points (for constructing the index-at-length and index-at-age) in the data-list used by TMB.

This function returns the data-list used by TMB
}
\description{
includeIntPoints
}
