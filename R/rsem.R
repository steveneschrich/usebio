
# Some notes:
# ideally:
# import_rsem(sample_table, formula, output_type=c("gene","tx2gene","tx"), formula=~1, gene_model = "")
#
# read files in


#' Import RSEM as ExpressionSet
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


#' Import RSEM as SummarizedExperiment
#'
#'  [SummarizedExperiment::SummarizedExperiment] object (Counts, Normalized Counts and TPM).
#' @param sample_table A data.frame of annotations and filenames. Can also be a path (default:".")
#' or a vector of filenames. If a data frame, should have one column (filenames) or have the variable
#' 'filename' for the filenames.
#' @param txIn logical - Transcript-level input (vs gene)
#' @param TxOut logical - Transcript-level output (vs gene)
#' @param formula A formula for DESeq (default is ~1)
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] representing the RSEM
#' data (with assays for counts and TPM).
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
import_rsem_as_DESeqDataSet <- function(sample_table=".", txIn = FALSE, txOut = FALSE, formula = ~1) {


  # sample_table could be a data.frame (of annotation), a list of files,
  # or a single path to files. Handle the various cases and end up with
  # a data frame with 'filename' column (possibly among others).
  if ( is.data.frame(sample_table)) {
    if ( !utils::hasName(sample_table, "filename") ) {
      if ( ncol(sample_table) > 1 ) {
        cli::cli_abort("Sample table must have a 'filename' column or only one column.")
      } else {
        colnames(sample_table) <- c("filename")
      }
    }
  } else {
    # the sample_table is just a list of files (or a directory)
    if ( length(sample_table == 1) && all(fs::is_dir(sample_table)) ) {
      sample_table <- tibble::tibble(
        filename = find_rsem_files(sample_table, which = ifelse(txIn, "isoforms","genes"))
      )
    } else {
      sample_table <- tibble::tibble(
        filename = sample_table
      )
    }
  }

  # We need a sample name, which can be inferred or already provided.
  if ( !utils::hasName(sample_table, "sample") )
    sample_table$sample <- infer_rsem_samplename(sample_table$filename)

  # After all the preprocessing, we can import the results using the
  # filenames and sample names.
  txi.rsem <- import_rsem(
    files = sample_table$filename, sample_names = sample_table$sample,
    txIn = txIn, txOut = txOut
  )

  # There are DESeq2-specific things to do (convert to DESeq2 object and normalize)
  dds <- DESeq2::DESeqDataSetFromTximport(txi.rsem, sample_table, formula)
  dds <- DESeq2::estimateSizeFactors(dds)

  # And we want to add the TPM values as an extra assay to the object.
  SummarizedExperiment::assays(dds)$TPM <- txi.rsem$abundance

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
#' @param files A character vector of filenames to load. If not specified or if
#'  the input is a directory name, all appropriate files are found in the specified
#'  directory (or current directory).
#' @param sample_names Optional character vector of sample names to use for files. If not
#'  specified, the sample name from the filename will be inferred.
#' @param which Select either 'genes' or 'isoforms' for corresponding RSEM outputs
#'
#' @return A list (see [tximport::tximport()]) of matrices containing RSEM data.
#' @export
#'
#' @examples
#' \dontrun{
#' import_rsem(c("foo.txt","bar.txt"))
#' }
import_rsem <- function(files = ".", sample_names = NULL, txIn = FALSE, txOut = FALSE) {
                        #which = c("genes","isoforms")) {

  # Select either gene-level or isoform-level input/output
  which_input <- ifelse(txIn, "isoforms", "genes")
  # If files is length 1 and a directory, list the files.
  if ( length(files) == 1 && all(fs::is_dir(files) )) {
    files <- find_rsem_files(dir = files, which = which_input)
  }
  # If there were no sample names, use the basename of the given files
  if ( is.null(sample_names) )
    sample_names <- infer_rsem_samplename(files)

  # If we are using txIn=TRUE (transcript-level) and output is gene level (txOut=FALSE),
  # assemble the transcript/gene mapping by pre-loading all data to find the mappings.
  # This is super-inefficient.
  if ( txIn && !txOut ) {
    tx2gene <- get_tx2gene_mapping(files)
    x <- tximport::tximport(files, type="rsem", txIn =txIn, txOut = txOut, tx2gene = tx2gene)
  } else {
    # Import expression data (gene or isoform level)
    x <- tximport::tximport(files, type = "rsem", txIn = txIn, txOut = txOut)
  }

  # Add sample names
  x <- add_colnames(x, sample_names)

  #
  #x <- remove_empty_rows_from_rsem(x)
  x
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
#' @description Find RSEM files ending with appropriate extensions for loading.
#'
#' @details Identify files ending with `.genes.results` or `.isoforms.results` as RSEM
#' files suitable for loading (with full
#' pathname).
#'
#' @param dir A directory to look in
#' @param which Select either 'genes' or 'isoforms' for corresponding RSEM outputs
#' @return A character vector of full paths to RSEM output files
#' @export
#'
#' @examples
find_rsem_files <- function(dir, which = c("genes","isoforms")) {
  which <- match.arg(which, c("genes","isoforms"))
  list.files(path = dir, pattern = sprintf("*.%s.results", which), full.names = TRUE)
}



#' Extract sample name from RSEM filename
#'
#' @param s A string (or vector of strings) representing filenames
#'
#' @return A vector (same length as s) consisting of sample names (stripped filenames).
#' @export
#'
#' @examples
infer_rsem_samplename <- function(s) {
  basename(s) |>
    stringr::str_remove("\\.(genes|isoforms)\\.results$")
}

#' Generate a RSEM-based sample table based on a directory
#'
#' This function infers a set of samples that have RSEM output based on
#' a target directory (`dir`) and a RSEM type (`genes` or `isoforms`). The
#' sample name is then inferred from these files, resulting in a data frame
#' with `files` and `names` columns.
#'
#' @param dir The directory to look in for RSEM files.
#' @param which What RSEM output type (genes, isoforms) to identify
#'
#' @return A data.frame consisting of two columns `files` and `names`.
#' @export
#'
build_rsem_sample_table <- function(dir, which = c("isoforms","genes")) {
  which <- match.arg(which)

  tibble::tibble(
    files = find_rsem_files(dir, which = which),
    names = infer_rsem_samplename(files)
  )
}

#' Add column names to matrices in the list x
#'
#' @param x A list of matrices and other elements (ignored)
#' @param colnames A list of column names to set for all matrices
#'
#' @return
#' @export
#'
#' @examples
add_colnames <- function(x, colnames) {
  stopifnot(is.list(x))
  purrr::map(x, function(.x) {
    if ( is.matrix(.x) )
      colnames(.x) <- colnames
    .x
  })
}

get_tx2gene_mapping <- function(files = ".") {
  if ( length(files) == 1 && all(fs::is_dir(files))) {
    files <- find_rsem_files(files, which = "isoforms")
  }

  # For all input files, read in the first two columns (transcript_id and gene_id).
  # The goal is to create a tx2gene data frame suitable for tximport. This is
  # just a unique data frame.
  purrr::map(files, \(f) {
    readr::read_tsv(file = f, col_select=c("transcript_id","gene_id"), col_types="cc")
  }) |>
    purrr::list_rbind() |>
    dplyr::distinct()
}
