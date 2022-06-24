#' Read a GCT-formatted file into an ExpressionSet object.
#'
#' @param file The name of the file (GCT format) to read in.
#'
#' @return An ExpressionSet representing the GCT data.
#' @export
#'
#' @examples
#' eset<-readGCT("/tmp/data.gct")
#'
readGCT<-function(file) {
  require(Biobase)

  stopifnot(typeof(file)=="character")

  # A GCT is a matrix (after line 2), but with two identifier columns
  d<-read.table(file=file,
                header=TRUE,
                row.names=1,
                sep="\t",
                skip=2,
                as.is=TRUE,
                check.names=FALSE,
                comment.char="",
                na.strings="",
                quote="",
                blank.lines.skip=TRUE
  )

  vals<-data.matrix(d[,-1])
  descriptions<-data.frame(Description=d[,1],
                           row.names=rownames(d),
                           stringsAsFactors = FALSE
  )
  res<-Biobase::ExpressionSet(assayData=vals,
                              featureData=AnnotatedDataFrame(descriptions))

  return(res)
}


#' Write an ExpressionSet to a GCT-formatted output file.
#'
#' @param eset The ExpressionSet to write.
#' @param file The file to write to.
#'
#' @return NULL
#' @export
#'
#' @examples
#' writeGCT(eset, file="/tmp/exprs.gct")
#'
writeGCT<-function(eset, file) {
  stopifnot(class(eset)=="ExpressionSet")

  # GCT requires two columns: a NAME and Description. They will both be the
  # featureNames from the ExpressionSet.
  output<-data.frame(NAME=featureNames(eset),
                     Description=featureNames(eset),
                     exprs(eset),
                     stringsAsFactors = FALSE,
                     check.names = FALSE)

  # Create the header separately
  hdr<-data.frame(a=c("#1.2",dim(eset)[1]),b=c("",dim(eset)[2]))

  # This bit writes the header first, then appends the data. The warnings are
  # suppressed for the second part, since it's a different size. The first one
  # may fail due to file not being created, so the suppress may be ok.
  write.table(hdr, file=file, quote=F, sep="\t", row.names=FALSE, col.names=FALSE)
  suppressWarnings(
    write.table(output, file=file, append=TRUE, quote=F, sep="\t",row.names=FALSE)
  )


}

