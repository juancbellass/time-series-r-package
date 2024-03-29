% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Similarity_indices.R
\name{Pearsons.Measure}
\alias{Pearsons.Measure}
\title{PEARSON MEASURE}
\usage{
Pearsons.Measure(s1, s2)
}
\arguments{
\item{s1}{a numeric vector containing the first time series.}

\item{s2}{a numeric vector containing the second time series.}
}
\value{
a positive real number between 0 and 1. Where 0 is the maximum similarity and 1 the maximum dissimilarity.
}
\description{
Computes the Pearson correlation between two numerical vectors and return the result of 1-abs(Person).
}
