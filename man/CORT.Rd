% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Similarity_indices.R
\name{CORT}
\alias{CORT}
\title{TEMPORAL CORRELATION}
\usage{
CORT(s1, s2)
}
\arguments{
\item{s1}{a numeric vector containing the first time series.}

\item{s2}{a numeric vector containing the second time series.}

\item{k}{a equal to 3.1 by default.}
}
\value{
CORT temporal correlation between S1 and S2.
}
\description{
Computes the Temporal Correlation of Chouakria and Nagabhushan between two same length time series S1 and S2.
}
\references{
for more details read  Chouakria Douzal 2003.
}
