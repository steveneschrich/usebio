#' Annotate object indexed by ensembl
#'
#' @description Annotate an object with gene identifiers that are ensembl
#'
#' @details
#'
#' @note This function currently only handles a SummarizedExperiment or child object.
#' @param x An object (ExpressionSet, SummarizedExperiment, etc)
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#' @examples
annotate_ensembl <- function(x) {
  ensembl_ids <- rownames(x)

  stopifnot(all(stringr::str_starts(ensembl_ids,"ENSG")))

  if ( any(stringr::str_ends(ensembl_ids, "\\.\\d+")) ) {
    ensembl_ids <- stringr::str_remove(ensembl_ids, "\\.\\d+$")
  }

  # Get annotation for ensembl ids from R database
  gene_annotation <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL")
  )

  # Since it is likely not a 1:1, select the lowest ENTREZID entry associated with ENSEMBL
  gene_annotation <- gene_annotation |>
    dplyr::group_by(.data$ENSEMBL) |>
    dplyr::arrange(.data$ENTREZID) |>
    dplyr::slice_head(n=1)

  # Need to make sure the mapping went correctly.
  stopifnot(all(ensembl_ids== gene_annotation$ENSEMBL))

  # Create rowData annotation
  SummarizedExperiment::rowData(x) <- S4Vectors::cbind(
    SummarizedExperiment::rowData(x),
    gene_annotation
  )

  x

}
