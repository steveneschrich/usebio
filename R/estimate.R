#' Run the ESTIMATE algorithm on an ExpressionSet, and return the result as
#' an ExpressionSet.
#'
#' @param eset An affy gene expression dataset.
#' @param gene_symbols A list of gene symbols (if not supplied in annotation for ExpressionSet).
#' @return An expression set consisting of estimate scores.
#' @export
#'
#' @examples
#' \dontrun{
#' estimate_scores<-runEstimate(eset)
#' }
runEstimate<-function(eset, platform, ...) {
  stopifnot(!is.ExpressionSet(eset))

  if (missing(platform)) platform<-"affymetrix"

  # Transform the eset to index by gene symbol
  eset_bygene<-affyutils::summarize_by_symbol(eset, by = "max")

  eset_common_genes<-affyutils::subset_by(eset_bygene, subset = estimate::common_genes$GeneSymbol)

  # Estimate is file-based
  estimate_exprs_in_file <- tempfile(pattern="estimate", fileext=".gct")
  writeGCT(eset_common_genes, estimate_exprs_in_file)

  estimate_scores_file <-tempfile(pattern="estimate-out", fileext=".gct")
  estimateScore(input.ds=estimate_exprs_in_file,
                output.ds=estimate_scores_file,
                platform=platform)

  results<-readGCT(estimate_scores_file)

  # At the moment, estimateScore converts to a data.frame inside of it, thereby
  # bashing sample names (if they are not valid). We just replace the sample
  # names with the initial ones.
  Biobase::sampleNames(results)<-Biobase::sampleNames(eset)

  results
}
