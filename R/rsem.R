

#' Import RSEM as SummarizedExperiment
#'
#'  [SummarizedExperiment::SummarizedExperiment] object (Counts, Normalized Counts and TPM).
#' @param sample_table A data.frame of annotations and filenames. Can also be a path (default:".")
#' or a vector of filenames. If a data frame, should have one column (files) or have the variable
#' 'files' for the filenames. There can also be a `names` variable which will be used for
#' sample names.
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
  txi.rsem <- import_rsem_files(
    files = sample_table$filename, names = sample_table$sample,
    txIn = txIn, txOut = txOut
  )

  # There are DESeq2-specific things to do (convert to DESeq2 object and normalize)
  dds <- DESeq2::DESeqDataSetFromTximport(txi.rsem, sample_table, formula)
  dds <- DESeq2::estimateSizeFactors(dds)

  # And we want to add the TPM values as an extra assay to the object.
  SummarizedExperiment::assays(dds)$TPM <- txi.rsem$abundance

  # Add row and col data.

  dds
}


#' Verify that tximport list is valid
#'
#' @param x A list from the [tximport::tximport()] function
#'
#' @return A logical indicating of the tximport is valid.
#' @export
#'
valid_tximport <- function(x) {
  # Make sure x is a tximport list.
  res <- utils::hasName(x, "abundance") &&
    utils::hasName(x, "counts") && utils::hasName(x,"length")

  # Make sure dimension names are consistent across list elements
  res <- res &&
    all(colnames(x[["abundance"]])==colnames(x[["counts"]])) &&
    all(colnames(x[["abundance"]])==colnames(x[["length"]])) &&
    all(rownames(x[["abundance"]])==rownames(x[["counts"]])) &&
    all(rownames(x[["abundance"]])==rownames(x[["length"]]))

  res
}


#' Import tximport list to DESeq2 object
#'
#' Convert a tximport list of elements from [tximport::tximport()] into a
#' DESeq2 object [DESeq2::DESeqDataSet()] for downstream processing.
#'
#' @param x The tximport list from [tximport::tximport()]
#'
#' @return
#' @export
#'
tximport_to_DESeq2 <- function(x) {
  stopifnot(valid_tximport(x))

  # There are DESeq2-specific things to do (convert to DESeq2 object and normalize)
  xdeseq <- DESeq2::DESeqDataSetFromTximport(x, colData = x[["sample_table"]], ~1)
  xdeseq <- DESeq2::estimateSizeFactors(xdeseq)

  # And we want to add the TPM values as an extra assay to the object.
  SummarizedExperiment::assays(xdeseq)$TPM <- x[["abundance"]]


  xdeseq
}

#' Title
#'
#' @param x
#' @param gtf_url
#' @param which
#'
#' @return
#' @export
#'
#' @examples
annotate_rsem_rows <- function(x, gtf_url=NULL, which = c("gene","transcript")) {
  which <- match.arg(which)

  # Retrieve gene/transcript-level annotation from the gtf
  g <- as.data.frame(gtf_annotation(gtf_url, feature.type = which))

  g <- dplyr::rename(
    g,
    gene_seqnames = seqnames,
    gene_start = start,
    gene_end = end,
    gene_width = width,
    gene_strand = strand
  )
  # Set annotation rownames to identifier (gene or transcript)
  rownames(g) <- g[[ifelse(which=="gene","gene_id","transcript_id")]]
  stopifnot(
    all(rownames(g) %in% rownames(x)),
    all(rownames(x) %in% rownames(g))
  )
  # Rearrange gene annotation relative to the experimental object
  g <- g[rownames(x),]

  # Then add the annotation as row annotation
  SummarizedExperiment::rowData(x) <- g

  x
}

#' Title
#'
#' @param x
#' @param phenoData
#' @param featureData
#'
#' @return
#' @export
#'
#' @examples
tximport_to_ExpressionSet <- function(x, phenoData = NULL, featureData = NULL) {

  stopifnot(valid_tximport(x))

  # Make sure phenoData (if provided) matches the colnames
  stopifnot((!is.null(phenoData) && all(colnames(x[["abundance"]]) %in% rownames(phenoData))))
  # Make sure the featureData (if provided) matches
  stopifnot((!is.null(featureData) && all(rownames(x[["abundance"]]) %in% rownames(featureData))))

  # There is an environment of assayData to store measurements
  ad <- new.env()
  ad$TPM <- x[["abundance"]]
  ad$counts <- x[["counts"]]
  ad$length <- x[["length"]]

  if ( !is.null(phenoData)) {
    phenoData <- Biobase::AnnotatedDataFrame(
      phenoData[colnames(x[["abundance"]]),]
    )
  }
  if ( !is.null(featureData)) {
    featureData <- Biobase::AnnotatedDataFrame(
      featureData[rownames(x[["abundance"]]),]
    )
  }
  es <- Biobase::ExpressionSet(
    assayData = ad,
    phenoData = phenoData,
    featureData = featureData,
    annotation = "tximport"
  )
  es
}


#' Title
#'
#' @param x
#' @param col_data
#' @param row_data
#'
#' @return
#' @export
#'
#' @examples
tximport_to_SummarizedExperiment <- function(x,col_data = NULL, row_data = NULL) {
  stopifnot(valid_tximport(x))

  stopifnot(
    all(!is.null(col_data) && colnames(x[["abundance"]] %in% rownames(col_data))),
    all(!is.null(row_data) && rownames(x[["abundance"]] %in% rownames(row_data)))
  )

  if ( !is.null(col_data)) {
    col_data <- Biobase::AnnotatedDataFrame(
      col_data[colnames(x[["abundance"]]),]
    )
  }
  if ( !is.null(row_data)) {
    row_data <- Biobase::AnnotatedDataFrame(
      row_data[rownames(x[["abundance"]]),]
    )
  }
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = x[["counts"]], abundance = x[["abundance"]], length = x[["length"]]),
    rowData = row_data,
    colData = col_data

  )
}

#' Title
#'
#' @param sample_table
#' @param which
#' @param tx2gene
#' @param gene_annotation
#' @param importer
#' @return
#' @export
#'
#' @examples
import_rsem <- function(
    sample_table,
    which = c("gene","transcript","tx2gene"),
    tx2gene = NULL,
    gene_annotation = NULL,
    importer = NULL
) {

  which <- match.arg(which)

  x <- usebio::tximport_rsem(sample_table, which = which, tx2gene = tx2gene, importer=importer)
  x <- usebio::tximport_to_DESeq2(x)
  x <- usebio::annotate_rsem_rows(
    x, gtf_url = gene_annotation,
    which = ifelse(which=="transcript","transcript","gene")
  )

  x
}
#' Import RSEM data using a sample table
#'
#' Import data from files using a sample table (with files and optionally names),
#' either at the transcript or gene level.
#'
#' @param sample_table
#' @param which
#' @param tx2gene
#' @param importer
#' @return
#' @export
#'
tximport_rsem <- function(sample_table, which=c("gene","transcript","tx2gene"),tx2gene=NULL,importer=NULL) {
  tximport_rsem_files(
    files = sample_table[["files"]],
    names = sample_table[["names"]],
    txIn = which %in% c("transcript","tx2gene"),
    txOut = which %in% c("transcript"),
    tx2gene_gtf = tx2gene,
    importer = importer
  )
}

#' Import RSEM data using tximport
#'
#' @description Import raw RSEM data into a list of quantifications.
#'
#' @details This function takes a sample table containing the `files` field,
#' and uses [tximport::tximport()] to load the files. The result is a list
#' containing several elements: abundance, counts, length. These are matrices
#' representing the RSEM data for a list of samples.
#'
#' @note This function is focused on gene-level RSEM data, not transcript-level.
#'
#' @param files A character vector of filenames to load. If not specified or if
#'  the input is a directory name, all appropriate files are found in the specified
#'  directory (or current directory).
#' @param names Optional character vector of sample names to use for files. If not
#'  specified, the sample name from the filename will be inferred.
#' @param txIn (logical) See [tximport::tximport()], is transcript-level input expected.
#' @param txOut (logical) See [tximport::tximport()], is transcript-level output required.
#' @param tx2gene_gtf (string) If txOut is FALSE and txIn is TRUE, then the data must be
#'  summarized from transcript-level to gene level. The tx2gene mapping can be inferred
#'  from the input file(s), or if a gtf annotation file is available, use this. This
#'  parameter is for a URL/file associated with a gtf for this conversion.
#' @param importer A function (see [tximport::tximport()]) for reading in the rsem output
#'  files. Note that if NULL, options include using [arrow::read_tsv_arrow()], [readr::read_tsv()].
#'
#' @return A list (see [tximport::tximport()]) of matrices containing RSEM data. Note that
#' a `sample_table` list element will be added (data frame) with `files` and `names` as
#' columns (based on input data).
#'
#' @export
#'
#'
tximport_rsem_files <- function(
    files = ".",
    names = NULL,
    txIn = FALSE,
    txOut = FALSE,
    tx2gene_gtf = NULL,
    importer = NULL
) {

  # Select either gene-level or isoform-level input/output
  which_input <- ifelse(txIn, "transcript", "gene")

  # If files is length 1 and a directory, list the files.
  if ( length(files) == 1 && all(fs::is_dir(files) )) {
    files <- find_rsem_files(dir = files, which = which_input)
  }
  # If there were no sample names, use the basename of the given files
  if ( is.null(names) )
    names <- infer_rsem_samplename(files)

  # Importer: If arrow is installed use this
  if ( is.null(importer) ) {
    if ( rlang::is_installed("arrow") )
      importer <- arrow::read_tsv_arrow
    else
      cli::cli_alert_info("Package arrow not installed, consider using it for efficiency.")
  }

  # Import the data as a list
  x <- tximport::tximport(files, type = "rsem", txIn = txIn, txOut = txIn,
                          importer = importer)
  # Add sample names
  x <- add_colnames(x, names)

  # Special case: if txIn and !txOut, we are summarizing to gene level via
  # tximport. However, this requires a mapping.
  # assemble the transcript/gene mapping by pre-loading all data to find the mappings.
  # This is super-inefficient.
  if ( txIn && !txOut ) {
    if ( is.null(tx2gene_gtf) )
      tx2gene_mapping <- infer_tx2gene_mapping(files)
    else {
      tx2gene_mapping <- gtf_tx2gene(tx2gene_gtf)
    }
    x <- tximport::summarizeToGene(x, tx2gene = tx2gene_mapping)
  }

  # Include the sample table in the list.
  x[["sample_table"]] <- tibble::tibble(files = files, names = names)

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
find_rsem_files <- function(dir, which = c("gene","transcript")) {
  which <- match.arg(which)
  list.files(
    path = dir,
    pattern = sprintf(
      "*.%s.results",
      ifelse(which=="gene","genes","isoforms")
    ),
    full.names = TRUE
  )

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
build_rsem_sample_table <- function(dir, which = c("gene","transcript")) {
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

#' Infer tx2gene mapping for RSEM
#'
#' A simple approach of reading the tx2gene mapping directly from the
#' first RSEM input file provided. Not a great strategy but sometimes
#' this is all that is available.
#'
#' @param files A list of RSEM transcript-level files
#'
#' @return A tx2gene mapping (two column data frame, first is tx, second is gene).
#' @export
#'
infer_tx2gene_mapping <- function(files = ".") {
  if ( length(files) == 1 && all(fs::is_dir(files))) {
    files <- find_rsem_files(files, which = "transcript")
  }

  # For all input files, read in the first two columns (transcript_id and gene_id).
  # The goal is to create a tx2gene data frame suitable for tximport. This is
  # just a unique data frame.
  purrr::map(files[1], \(f) {
    readr::read_tsv(file = f, col_select=c("transcript_id","gene_id"), col_types="cc")
  }) |>
    purrr::list_rbind() |>
    dplyr::distinct()
}
