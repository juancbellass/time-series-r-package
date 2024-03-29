% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Similarity_indices.R
\name{DTW.Measure}
\alias{DTW.Measure}
\title{DYNAMIC TIME WARPING
Computes the Dynamic Time Warping distance between a pair of numeric time series using the 'DTWDistance' function from the 'TSdist' package.
We by setting the 'DTWDistance' function by using a  Sakoe-Chiba windowing constraint with a window equal to the 10% of the time series length.}
\usage{
DTW.Measure(s1, s2)
}
\arguments{
\item{s1}{a numeric vector containing the first time series.}

\item{s2}{a numeric vector containing the second time series.}
}
\value{
the computed distance between the pair of series.
}
\description{
DYNAMIC TIME WARPING
Computes the Dynamic Time Warping distance between a pair of numeric time series using the 'DTWDistance' function from the 'TSdist' package.
We by setting the 'DTWDistance' function by using a  Sakoe-Chiba windowing constraint with a window equal to the 10% of the time series length.
}
\references{
for more details read the 'TSdist' package documentation.
}
