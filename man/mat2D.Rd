% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mat2D.R
\name{mat2D}
\alias{mat2D}
\title{Generates Table of Mean 1D CNV Profile}
\usage{
mat2D(
  cncf,
  bin.size = 10,
  amp.tcn = 4,
  clustering.prep = TRUE,
  progress.bar = TRUE
)
}
\arguments{
\item{cncf}{Input cncf file}

\item{bin.size}{Size of bin in mb}

\item{amp.tcn}{Lowest tcn for AMP}

\item{clustering.prep}{Generate chromosome level features for clustering}

\item{progress.bar}{Display progress bar}
}
\value{
List of components \code{bins} (chromosome, start, and end of bins),
        and \code{matrix} (sample vs. bin matrix of copy numbers).
}
\description{
Generates Table of Mean 1D CNV Profile
}
\examples{
z <- mat2D(cncf = system.file('extdata','example_with_ploidy.cncf', package = 'facetsHeatmap'))
names(z)
head(z$matrix)

}
