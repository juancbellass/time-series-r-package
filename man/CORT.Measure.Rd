% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Similarity_indices.R
\name{CORT.Measure}
\alias{CORT.Measure}
\title{TEMPORAL CORRELATION MEASURE}
\usage{
CORT.Measure(s1, s2)
}
\arguments{
\item{s1}{a numeric vector containing the first time series.}

\item{s2}{a numeric vector containing the second time series.}

\item{k}{a equal to 3.1 by default.}
}
\value{
a positive real number between 0 and 1. Where 0 is the maximum similarity and 1 the maximum dissimilarity.
}
\description{
Computes the Temporal Correlation of Chouakria and Nagabhushan between two same length time series s1 and s2, and return the result of 1-abs(CORT).
}
\references{
for more details read  Chouakria Douzal 2003.
}
