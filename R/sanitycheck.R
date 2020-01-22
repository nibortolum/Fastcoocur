#' Check that the input matrix meet conditions to compute cooccurences
#'
#' @param spp_site_mat The species matrix (apecies (rows) x sites (column))
#' @param strict TRUE or FALSE. If TRUE, it will throw an error at the first encountered problem and stop code execution.
#' If FALSE, if will perform all the tests and throw warnings.
#'
#' @return A list of warnings.
#' @export
#'
#' @examples
sanity_check <- function(spp_site_mat, strict = FALSE){

  all.good <- TRUE
  if(!is.matrix(spp_site_mat)) {
    err <- paste("spp_site_mat must be a matrix. Try coercing with as.matrix()\n")
    if(strict) stop(err)
    warning(err)
    all.good <- FALSE
  }

  if(is.null(row.names(spp_site_mat))) {
    warning("spp_site_mat as no rownames. fqst_cooccur will generate dummy names\n")
  }

  if(any(rowSums(spp_site_mat) == 0)){
    err <- paste(" somes species are not present in any samples. cooccurence computation is hence impossible. Please remove them.\n")
    if(strict) stop(err)
    warning(err)
    all.good <- FALSE
  }

  if(all.good) cat("Tests passed successfully\n")

}
