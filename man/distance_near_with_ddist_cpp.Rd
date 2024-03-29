% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{distance_near_with_ddist_cpp}
\alias{distance_near_with_ddist_cpp}
\title{Calculate thresheld pairwise distance}
\usage{
distance_near_with_ddist_cpp(locs, locs_grid, radius, n_neighbors = 68L)
}
\arguments{
\item{locs}{An N x 2 matrix of spatial locations}

\item{locs_grid}{An N_grid x 2 matrix of spatial grids}

\item{radius}{The thresholded radius}

\item{n_neighbors}{The expected number of neighbors based on the MRA grid}
}
\value{
The thresheld pairwise distance
}
\description{
Calculate thresheld pairwise distance
}
