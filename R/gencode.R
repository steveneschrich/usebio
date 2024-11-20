#' gtf
#'
#' @description Routines for interacting with gtf gene annotation files
#'
#' @details GTF is a gene annotation format accessible from EBI/Ensembl. These
#' routines are wrappers around the excellent [rtracklayer] package. The goal
#' is merely to make working with these files a bit easier, in the mode of
#' the [usethis] package.
#'
#' The GTF files (for gencode) have information at the gene, transcript and
#' exon levels. It is possible to get an object with all of these, or with
#' only data at one of these levels. These functions make the process slightly
#' easier.
#'
#' *Caching*: Since the files in question are large, one of the things that
#' this package does is provide caching (via BiocFileCache) of the files for
#' subsequent loads. Hopefully, it is often the case that the same file would
#' be needed multiple times.
#'
#' Finally, although not very well thought out, there is a function (import_gencode_gtf)
#' which attempts to impose some logic on the process of inferring the right gencode
#' file. It is very much a work in progress.
#'
#'
#' @name gtf_gencode
NULL

#' Import a gencode gtf
#'
#' @param version A version number (assuming human)
#' @param url A URL to use for direct file access. Note all other parameters
#'  are ignored if this is set.
#' @param species The species (currently human) to import
#' @param gencode_base_url The base URL of the genocode download (EBI)
#'
#' @return A GRanges object representing the genocode gtf (see [rtracklayer::import.gff()]).
#' @export
#'
import_gencode_gtf <- function(
    url = NULL,
    version = "v32",
    species = c("human"),
    gencode_base_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode"
) {

  species <- match.arg(species)

  # Attempt to normalize version information
  if ( is.numeric(version) ) {
    numeric_version <- version
  } else {
    numeric_version <- regmatches(version, regexpr("[0-9]+$",version))
    numeric_version <- as.numeric(numeric_version)
  }

  # If we don't have a url, build it here.
  if ( is.null(url) ) {
    gencode_dir <- sprintf(
      "%s/%s/release_%d",
      gencode_base_url,
      ifelse(species=="human", "Gencode_human", "Gencode_human"),
      numeric_version
    )
    url <- sprintf("%s/gencode.v%d.annotation.gtf.gz",gencode_dir, numeric_version)
  }


  cli::cli_progress_step("Importing gtf to GRanges format")
  gtf_annotation(url)

}


#' Create tx2gene mapping from gtf
#'
#' @param gtf A gtf object to extract mapping
#'
#' @return Data frame with two columns (transcript_id and gene_id), suitable
#'  for use in [tximport::tximport()].
#' @export
#'
tx2gene <- function(url, ...) {

  t2g <- gtf_transcript_annotation(url, ...)

  dplyr::distinct(as_tibble.DFrame(t2g), transcript_id, gene_id)
}

#' Extract gene-level annotation from gtf
#'
#' @param gtf
#'
#' @return
#' @export
#'
#' @examples
gtf_gene_annotation <- function(url, ...) {
  gtf_annotation(url, feature.type="gene")
}

#' Extract transcript-level annotation from gtf
#'
#' @param gtf
#'
#' @return
#' @export
#'
#' @examples
gtf_transcript_annotation <- function(url, ...) {
  gtf_annotation(url, feature.type="transcript")
}

#' Title
#'
#' @param url
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
gtf_exon_annotation <- function(url, ...) {
  gtf_annotation(url, feature.type="exon")
}

#' Title
#'
#' @param url
#' @param feature.type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
gtf_annotation <- function(url, feature.type = NULL, ...) {
  g <- get_gtf(url)
  rtracklayer::import(g, format="gtf", feature.type = feature.type, ...)
}
#' Title
#'
#' @param url
#'
#' @return
#' @export
#'
#' @examples
get_gtf <- function(url) {
  # We utilize Bioconductor caching for the gtf
  bfc <- BiocFileCache::BiocFileCache(
    cache = rappdirs::user_cache_dir(appname="org.R-project.R/R/usebio"),
    ask = FALSE
  )

  cli::cli_progress_step("Retrieving gtf from {url}")
  BiocFileCache::bfcrpath(bfc, url)

}


# Some notes
#AnnotationHub::query(ah, c("TxDb","Gencode","v32", "hg38","gtf"))
#gtf <- ah[["AH75191"]]
# NB: This is a TxDb object
# Could also just have the range from rtracklayer::import. There is
# bioconductor caching of some type.



