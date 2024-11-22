# RNAseq

#' Title
#'
#' @param x
#' @param normalized
#'
#' @return
#' @export
#'
#' @examples
lnc <- function(x, normalized = FALSE) {
  stopifnot(is(x, "DESeqDataSet"))
  log2(DESeq2::counts(x, normalized=normalized)+0.1)
}
