

#' Install the ESTIMATE software (it is in a custom location)
#'
#' @return Invisible NULL (from install.packages).
#' @export
#'
#' @examples
install_estimate<-function() {

  rforge <- "http://r-forge.r-project.org"
  utils::install.packages("estimate", repos=rforge, dependencies=TRUE)
}




#' Install the XCell software (it is in a custom location)
#'
#' @return Invisible NULL (from install.packages).
#' @export
#'
#' @examples
install_xcell<-function() {

  devtools::install_github('dviraran/xCell')
}
