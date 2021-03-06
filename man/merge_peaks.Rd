% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Merge_Peaks.R
\name{merge_peaks}
\alias{merge_peaks}
\title{Given a data frame of peaks, merge peaks that are close together.}
\usage{
merge_peaks(x, min.distance = 5)
}
\arguments{
\item{x}{A data frame of peaks wich columns `Peak`, `Start`, and `End`.}

\item{min.distance}{The minimum distance allowed between two peaks.}
}
\value{
A dataframe with as many or fewer peaks.
}
\description{
Given a data frame of peaks, merge peaks that are close together.
}
