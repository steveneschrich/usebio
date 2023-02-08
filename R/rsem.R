
#' Title
#'
#'  [Biobase::ExpressionSet] object (Counts, Normalized Counts and TPM).
#' @param sample_table A data.frame of annotations and filenames
#' @param formula A formula for DESeq (default is ~1)
#' @param use_log Should the data be log2 transformed?
#'
#' @return A [Biobase::ExpressionSet] representing the RSEM data (with assayData for
#' exprs (normalized counts), counts (raw counts) and TPM).
#'
#' @export
#'
#' @examples
import_rsem_as_ExpressionSet <- function(sample_table, formula = ~1, use_log = TRUE) {
  stopifnot(utils::hasName(sample_table, "filename"))

  txi.rsem <- import_rsem(sample_table$filename)
  dds <- normalize_rsem(txi.rsem, sample_table, formula)

  # NB: The RSEM data does not have column names, so we need to label these
  # and use the sample_table sample names (or basename) as the rownames.
  if ( !utils::hasName(sample_table, "sample")) {
    sample_table <- sample_table |>
        dplyr::mutate(sample = stringr::str_remove(basename(sample_table$filename), ".genes.results"))
  }

  ad <- new.env()
  ad$TPM <- txi.rsem$abundance |> magrittr::set_colnames(sample_table$sample)
  ad$exprs <- DESeq2::counts(dds, normalized = TRUE) |> magrittr::set_colnames(sample_table$sample)
  ad$counts <- DESeq2::counts(dds, normalized = FALSE) |> magrittr::set_colnames(sample_table$sample)

  if ( use_log ) {
    ad$TPM <- log2(ad$TPM + 1)
    ad$exprs <- log2(ad$exprs + 1)
    ad$counts <- log2(ad$counts + 1)
  }
  es <- Biobase::ExpressionSet(
    assayData = ad,
    phenoData = Biobase::AnnotatedDataFrame(
      sample_table |>
        tibble::column_to_rownames("sample") |>
        as.data.frame()
      ),
    #    featureData = ?,
    annotation = "gene"
  )



  es
}


#' Title
#'
#'  [SummarizedExperiment::SummarizedExperiment] object (Counts, Normalized Counts and TPM).
#' @param sample_table A data.frame of annotations and filenames
#' @param formula A formula for DESeq (default is ~1)
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] representing the RSEM
#' data (with assays for counts and TPM).
#'
#' @export
#'
#' @examples
import_rsem_as_DESeqDataSet <- function(sample_table, formula = ~1) {
  stopifnot(utils::hasName(sample_table, "filename"))

  txi.rsem <- import_rsem(sample_table$filename)
  dds <- normalize_rsem(txi.rsem, sample_table, formula)

  # NB: The RSEM data does not have column names, so we need to label these
  # and use the sample_table sample names (or basename) as the rownames.
  if ( !utils::hasName(sample_table, "sample")) {
    sample_table <- sample_table |>
      dplyr::mutate(sample = stringr::str_remove(basename(sample_table$filename), ".genes.results"))
  }

  SummarizedExperiment::assays(dds)$TPM <- txi.rsem$abundance  |>
    magrittr::set_colnames(sample_table$sample)



  dds

}



#' Import RSEM data using tximport
#'
#' @description Import raw RSEM data into a list of quantifications.
#'
#' @details This function takes a sample table containing the `filename` field,
#' and uses [tximport::tximport()] to load the files. The result is a list
#' containing several elements: abundance, counts, length. These are matrices
#' representing the RSEM data for a list of samples.
#'
#' @note This function is focused on gene-level RSEM data, not transcript-level.
#'
#' @param files A character vector of filenames to load.
#'
#' @return A list (see [tximport::tximport()]) of matrices containing RSEM data.
#' @export
#'
#' @examples
#' \dontrun{
#' import_rsem(c("foo.txt","bar.txt"))
#' }
import_rsem <- function(files) {

 tximport::tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE) |>
    remove_empty_rows_from_rsem()
}

#' Title
#'
#' @param rsem RSEM object (list with counts, etc)
#' @param sample_table Sample table (metadata) for samples.
#' @param formula Formaula for DESeq (default ~1)
#'
#' @return
#' @export
#'
#' @examples
normalize_rsem <- function(rsem, sample_table, formula = ~1) {
  DESeq2::DESeqDataSetFromTximport(rsem, sample_table, formula) |>
    DESeq2::estimateSizeFactors()

}
#' Removes genes with no observations from RSEM object
#'
#' @description Removes rows (genes) that have length of 0 for any of the samples.
#'
#' @details There appears to be a need in DESeq2 to ensure that no gene have an estimated
#' gene length of 0. This routine filters out those rows from an RSEM object.
#'
#' @param x An RSEM object
#' @param remove.lengths (TRUE) Remove empty lengths
#' @param min.zero.obs (Inf) Number of zero observations to trigger removing row
#' @return An RSEM object with empty rows removed.
#' @export
#'
#' @examples
remove_empty_rows_from_rsem <- function(x, remove.lengths = TRUE, min.zero.obs = Inf) {
  haszero <- apply(x$length,1,min)==0
  too_many_zero_obs <- rowSums(x$counts==0) > min.zero.obs

  rows_to_remove <- haszero | (is.finite(min.zero.obs) & too_many_zero_obs)

  res<-x
  res$length <- x$length[!rows_to_remove,]
  res$abundance <- x$abundance[!rows_to_remove,]
  res$counts <- x$counts[!rows_to_remove,]

  res
}
#' Find RSEM files to load
#'
#' @description Identify files ending with `.genes.results` as RSEM files suitable for loading (with full
#' pathname).
#'
#' @param dir A directory to look in
#'
#' @return
#' @export
#'
#' @examples
find_rsem_files <- function(dir) {
  list.files(path = dir, pattern = "*.genes.results", full.names = TRUE)
}
