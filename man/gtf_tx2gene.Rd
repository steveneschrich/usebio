% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gencode.R
\name{gtf_tx2gene}
\alias{gtf_tx2gene}
\title{Create tx2gene mapping from gtf}
\usage{
gtf_tx2gene(url, ...)
}
\arguments{
\item{gtf}{A gtf object to extract mapping}
}
\value{
Data frame with two columns (transcript_id and gene_id), suitable
for use in \code{\link[tximport:tximport]{tximport::tximport()}}.
}
\description{
Create tx2gene mapping from gtf
}
