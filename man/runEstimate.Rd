% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{runEstimate}
\alias{runEstimate}
\title{Run the ESTIMATE algorithm on an ExpressionSet, and return the result as
an ExpressionSet.}
\usage{
runEstimate(eset, platform, ...)
}
\arguments{
\item{eset}{An affy gene expression dataset.}

\item{platform}{If missing, affymetrix}

\item{...}{Any other parameters to pass to estimate}
}
\value{
An expression set consisting of estimate scores.
}
\description{
Run the ESTIMATE algorithm on an ExpressionSet, and return the result as
an ExpressionSet.
}
\examples{
\dontrun{
estimate_scores<-runEstimate(eset)
}
}
