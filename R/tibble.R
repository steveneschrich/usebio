

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
#' @describeIn as_DFrame.tibble Use a data frame
#' @export
as_DFrame.data.frame <- function(x) {
  as_DFrame.tibble(x)
}

#' @export
as_DFrame <- function(x) {
  UseMethod("as_DFrame")
}


#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble


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
as_tibble.DFrame <- function(x) {
  # Strip out rownames before conversion
  rn <- rownames(x)
  if ( is.null(rn) ) rn <- ""
  rownames(x) <- NULL
  x <- cbind(.rownames = rn, as.data.frame(x))
  tibble::as_tibble(x)
}

#' Convert GRanges object to tibble
#'
#' @param x A [GenomicRanges::GRanges()] object
#'
#' @return A tibble representing `x`
#' @export
#'
as_tibble.GRanges <- function(x) {
  tibble::as_tibble(as.data.frame(x))
}

#' Extract column metadata as a tibble
#'
#' @description Extract the column metadata from an object and return
#' the result as a tibble.
#'
#' @param x The object to extract column metadata from
#'
#' @return A tibble containing column metadata from `x`.
#' @export
#' @docType methods
#' @rdname col_data-methods
#'
setGeneric("col_data", function(x) standardGeneric("col_data"))

#' Set column metadata as a tibble
#'
#' @description Some objects support column-oriented metadata. This function
#' will store this metadata with the object.
#'
#' @param x The tibble of metadata to store.
#'
#' @return The object with column metadata included.
#' @export
#' @docType methods
#' @rdname col_data-set
setGeneric("col_data<-", function(x, value) standardGeneric("col_data<-"))

#' @rdname col_data-methods
#' @aliases col_data,SummarizedExperiment,ANY-method
setMethod("col_data", "SummarizedExperiment",  function(x) {
  as_tibble(SummarizedExperiment::colData(x))
})

#' @rdname col_data-set
#' @aliases `col_data<-`,SummarizedExperiment,ANY-method
setMethod("col_data<-", "SummarizedExperiment", function(x, value) {
  set_col_data(x, value)
})

#' Title
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
set_col_data <- function(x, value) {
  SummarizedExperiment::colData(x) <- as_DFrame.tibble(value)
  x
}

#' Extract row data from an object
#'
#' @param x An object to extract row data from
#' @return A tibble of row data
#' @docType methods
#' @rdname row_data
#' @export
setGeneric("row_data", function(x) standardGeneric("row_data"))

#' Set row data on an object
#'
#' @description Some objects support row-oriented metadata. This function
#' will store this metadata with the object.
#'
#' @param x The object to set row data on
#' @param value The row data
#' @return The object with row data set
#' @export
#' @docType methods
#' @rdname row_data-set
setGeneric("row_data<-", function(x, value) standardGeneric("row_data<-"))

#' @rdname row_data
#' @aliases row_data,SummarizedExperiment,ANY-method
setMethod("row_data", "SummarizedExperiment", function(x) {
  as_tibble(SummarizedExperiment::rowData(x))
})

#' @rdname row_data-set
#' @aliases `row_data<-`,SummarizedExperiment,ANY-method
setMethod("row_data<-", "SummarizedExperiment", function(x, value) {
  set_row_data(x, value)
})

#' Title
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
set_row_data <- function(x, value) {
  SummarizedExperiment::rowData(x) <- as_DFrame.tibble(value)
  x
}
#' Title
#'
#' @param x
#' @param coldata
#'
#' @return
#' @export
#'
#' @examples
add_col_data <- function(x, coldata) {
  tdf <- usebio::col_data(x)
  if ( !utils::hasName(coldata, ".rownames") ) {
    stopifnot(!is.null(rownames(coldata)))
    coldata <- tibble::rownames_to_column(coldata, ".rownames")
  }
  tdf <- dplyr::left_join(tdf, coldata, by=".rownames",relationship="one-to-one")

  col_data(x) <- tdf
  x
}
#' Title
#'
#' @param x
#' @param rowdata
#'
#' @return
#' @export
#'
#' @examples
add_row_data <- function(x, rowdata) {
  tdf <- usebio::row_data(x)
  if ( !utils::hasName(rowdata, ".rownames") ) {
    stopifnot(!is.null(rownames(rowdata)))
    rowdata <- tibble::rownames_to_column(rowdata, ".rownames")
  }
  tdf <- dplyr::left_join(tdf, rowdata, by=".rownames",relationship="one-to-one")

  row_data(x) <- tdf
  x
}

