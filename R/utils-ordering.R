#' sort_idents
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the helpers family.
#' @param vec (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
sort_idents <- function(vec) {
  suffixes <<- sub(".*_", "", vec)
  
  # Determine if suffixes are numeric or alphabetic
  if (all(grepl("^[0-9]+$", suffixes))) {
    sorted_vec <<- vec[order(as.numeric(suffixes))]
  } else {
    sorted_vec <<- vec[order(suffixes)]
  }
  
  return(factor(rev(sorted_vec), levels = rev(sorted_vec)))
}


#' sort_by_reference
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the helpers family.
#' @param vec (auto) parameter
#' @param ref (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
sort_by_reference <- function(vec, ref) {
  # Find the intersection of vec and ref
  common_elements <- intersect(ref, vec)
  
  # Order vec according to the position of common elements in ref
  sorted_vec <- vec[order(match(vec, common_elements, nomatch = Inf))]
  
  # Remove elements that are not in common_elements
  sorted_vec <- sorted_vec[sorted_vec %in% common_elements]

  return(factor(rev(sorted_vec), levels = rev(sorted_vec)))
}
