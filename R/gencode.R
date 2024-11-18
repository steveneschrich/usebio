#' gtf
#'
#' @description Routines for interacting with gtf gene annotation files
#'
#' @details GTF is a gene annotation format accessible from EBI/Ensembl
#'
#' @name gtf
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

  # We utilize Bioconductor caching for the gtf
  bfc <- BiocFileCache::BiocFileCache(
    cache = rappdirs::user_cache_dir(appname="org.R-project.R/R/usebio"),
    ask = FALSE
  )

  cli::cli_progress_step("Retrieving gtf from {url}")
  gtf_cache <- BiocFileCache::bfcrpath(bfc, url)
  cli::cli_progress_step("Importing gtf to GRanges format")
  gtf_result <- rtracklayer::import(gtf_cache, format="gtf")


  gtf_result
}


#' Create tx2gene mapping from gtf
#'
#' @param gtf A gtf object to extract mapping
#'
#' @return Data frame with two columns (transcript_id and gene_id), suitable
#'  for use in [tximport::tximport()].
#' @export
#'
tx2gene <- function(gtf) {
  x <- as_tibble(gtf)

  dplyr::distinct(
    gtf_tx_annotation(x),
    transcript_id,
    gene_id
  )
}

#' Extract gene-level annotation from gtf
#'
#' @param gtf
#'
#' @return
#' @export
#'
#' @examples
gtf_gene_annotation <- function(gtf) {
  x <- as_tibble(grf)

 dplyr::filter(x, type == "gene")
}

#' Extract transcript-level annotation from gtf
#'
#' @param gtf
#'
#' @return
#' @export
#'
#' @examples
gtf_tx_annotation <- function(gtf) {
  x <- as_tibble(gtf)

  dplyr::filter(x, type == "transcript")
}

# Some notes
#AnnotationHub::query(ah, c("TxDb","Gencode","v32", "hg38","gtf"))
#gtf <- ah[["AH75191"]]
# NB: This is a TxDb object
# Could also just have the range from rtracklayer::import. There is
# bioconductor caching of some type.



