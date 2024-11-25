% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tibble.R
\docType{methods}
\name{row_data<-}
\alias{row_data<-}
\alias{row_data<-,SummarizedExperiment-method}
\alias{`row_data<-`,SummarizedExperiment,ANY-method}
\title{Set row data on an object}
\usage{
row_data(x) <- value

\S4method{row_data}{SummarizedExperiment}(x) <- value
}
\arguments{
\item{x}{The object to set row data on}

\item{value}{The row data}
}
\value{
The object with row data set
}
\description{
Some objects support row-oriented metadata. This function
will store this metadata with the object.
}
