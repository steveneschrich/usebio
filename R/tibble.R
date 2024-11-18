

#' Convert a tibble to a DFrame
#'
#' @details The tibble and DFrame are similar, but there are important
#' differences. This function does not handle all of them (yet). The
#' biggest thing it handles is managing the conversion of rownames,
#' which can be non-unique in DFrame but not in tibble/data.frame. It
#' is assumed (if using [as_tibble.DFrame()]) that the rownames are stored
#' in the data frame calumn `.rownames`. If it is not present, rownames are
#' not set. If non-unique, DFrame can handle this case by default.
#'
#' @param x A tibble (or data.frame)
#'
#' @return A [S4Vectors::DataFrame()] object coerced from `x`.
#' @export
#'
as_DFrame.tibble <- function(x) {
  rn <- x[[".rownames"]]
  if (all(rn=="")) rn <- NULL
  x[[".rownames"]]<-NULL
  S4Vectors::DataFrame(x, row.names=rn, check.names=FALSE)
}

#' Convert a DFrame to tibble
#'
#' @details There are a few differences between a tibble and DFrame, some
#'  of which are more involved. This function does not (yet) handle the
#'  complex cases.
#'
#'  rownames: tibbles don't like rownames (or dplyr), so we store the rownames
#'   in a variable called `.rownames` in the data frame. Should there be
#'   a need, see [tibble::column_to_rownames()] or [as_DFrame.tibble()] for
#'   options of getting the rownames back.
#' @param x A [S4Vectors::DFrame()] object
#'
#' @return A tibble representing `x`
#' @export
#'
as_tibble.DFrame <- function(x) {
  # Strip out rownames before conversion
  rn <- rownames(x)
  if ( is.null(rn) ) rn <- ""
  rownames(x) <- NULL
  cbind(.rownames = rn, as.data.frame(x))
}


setGeneric("col_data", function(x) standardGeneric("col_data"))
setGeneric("col_data<-", function(x, value) standardGeneric("col_data<-"))
setMethod("col_data", "SummarizedExperiment",  function(x) {
  as_tibble(SummarizedExperiment::colData(x))
})
setMethod("col_data<-", "SummarizedExperiment", function(x, value) {
  SummarizedExperiment::colData(x) <- as_DFrame.tibble(value)
  x
})

setGeneric("as_DFrame", function(x) standardGeneric("as_DFrame"))
setMethod("as_DFrame", "tbl_df", as_DFrame.tibble)
setGeneric("row_data", function(x) standardGeneric("row_data"))
setGeneric("row_data<-", function(x, value) standardGeneric("row_data<-"))
setMethod("row_data", "SummarizedExperiment", function(x) {
  as_tibble(SummarizedExperiment::rowData(x))
})
setMethod("row_data<-", "SummarizedExperiment", function(x, value) {
  SummarizedExperiment::rowData(x) <- as_DFrame.tibble(value)
  x
})
setGeneric("as_tibble", function(x) standardGeneric("as_tibble"))
setMethod("as_tibble","DFrame", as_tibble.DFrame)
setMethod("as_tibble", "GRanges", as_tibble.DFrame)

