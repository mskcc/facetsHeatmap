% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/showHeatmap.R
\name{showHeatmap}
\alias{showHeatmap}
\title{Generate heatmap of Copy-Number Profiles}
\usage{
showHeatmap(
  z,
  clusterZ = TRUE,
  K = 2,
  reorder.hc = TRUE,
  tcn.max = 6,
  relative = FALSE,
  tcn.min = -4,
  margin = c(0.15, 0.4, 0.5, 1.2),
  cex = 0.6,
  legend.title = "tcn",
  useLOH = FALSE,
  feature_output = NULL,
  show.sid = FALSE,
  sid.size = 2,
  ...
)
}
\arguments{
\item{z}{List of components \code{bins}, \code{matrix}, \code{dat.tcn},
and \code{dat.lcn}}

\item{clusterZ}{Cluster samples using abbreviated CNA features}

\item{K}{Number of clusters}

\item{reorder.hc}{Reorder samples by hierarchical clustering after grouping by
cluster}

\item{tcn.max}{Maximum total copy number to display}

\item{relative}{Use relative copy numbers}

\item{tcn.min}{If \code{Relative}, lower bound of relative tcn loss to show}

\item{margin}{Parameter to \code{theme(margin = ...)} for heatmap}

\item{legend.title}{Tile of legend}

\item{useLOH}{Show LOH heatmap instead of tcn}

\item{feature_output}{File name for sample-feature matrix; if \code{NULL},
matrix is not saved}

\item{show.sid}{Display sample IDs}

\item{sid.size}{Sample ID size}

\item{...}{Other parameters to \code{oneD(...)}}
}
\description{
Generate heatmap of Copy-Number Profiles
}
