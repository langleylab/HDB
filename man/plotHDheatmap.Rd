% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDB.R
\name{plotHDheatmap}
\alias{plotHDheatmap}
\title{Plot HDB sigmas (heatmap)}
\usage{
plotHDheatmap(hdb, group, rng = NULL)
}
\arguments{
\item{hdb}{the result of \code{HDB()}}

\item{group}{character containing the name of the \code{HDB} results list}

\item{rng}{vector containing a range in the form \code{c(a, b)} where
\code{a} and \code{b} contain the lower and upper bound of the scale.
Used to compare different heatmaps using the same scale.}
}
\value{
a heatmap with pairwise sigma values for the HDB calculation within
   \code{group}.
}
\description{
Plots the sigma values for every HDB comparison in a heatmap
}
