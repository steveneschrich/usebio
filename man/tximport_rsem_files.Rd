% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsem.R
\name{tximport_rsem_files}
\alias{tximport_rsem_files}
\title{Import RSEM data using tximport}
\usage{
tximport_rsem_files(
  files = ".",
  names = NULL,
  txIn = FALSE,
  txOut = FALSE,
  tx2gene_gtf = NULL,
  importer = NULL
)
}
\arguments{
\item{files}{A character vector of filenames to load. If not specified or if
the input is a directory name, all appropriate files are found in the specified
directory (or current directory).}

\item{names}{Optional character vector of sample names to use for files. If not
specified, the sample name from the filename will be inferred.}

\item{txIn}{(logical) See \code{\link[tximport:tximport]{tximport::tximport()}}, is transcript-level input expected.}

\item{txOut}{(logical) See \code{\link[tximport:tximport]{tximport::tximport()}}, is transcript-level output required.}

\item{tx2gene_gtf}{(string) If txOut is FALSE and txIn is TRUE, then the data must be
summarized from transcript-level to gene level. The tx2gene mapping can be inferred
from the input file(s), or if a gtf annotation file is available, use this. This
parameter is for a URL/file associated with a gtf for this conversion.}

\item{importer}{A function (see \code{\link[tximport:tximport]{tximport::tximport()}}) for reading in the rsem output
files. Note that if NULL, options include using \code{\link[arrow:read_delim_arrow]{arrow::read_tsv_arrow()}}, \code{\link[readr:read_delim]{readr::read_tsv()}}.}
}
\value{
A list (see \code{\link[tximport:tximport]{tximport::tximport()}}) of matrices containing RSEM data. Note that
a \code{sample_table} list element will be added (data frame) with \code{files} and \code{names} as
columns (based on input data).
}
\description{
Import raw RSEM data into a list of quantifications.
}
\details{
This function takes a sample table containing the \code{files} field,
and uses \code{\link[tximport:tximport]{tximport::tximport()}} to load the files. The result is a list
containing several elements: abundance, counts, length. These are matrices
representing the RSEM data for a list of samples.
}
\note{
This function is focused on gene-level RSEM data, not transcript-level.
}
