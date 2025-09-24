#' make_folder
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the helpers family.
#' @param folder_path (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
make_folder <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    # Create the folder
    dir.create(folder_path)
  }
}


#' strip_if_contains
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the helpers family.
#' @param vector (auto) parameter
#' @param search (auto) parameter
#' @param strip (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
strip_if_contains <- function(vector, search, strip) {
  # Initialize an empty vector to store results
  result_vector <- vector
  
  # Loop through each element of the vector
  for (i in seq_along(vector)) {
    if (grepl(search, vector[i])) {
      result_vector[i] <- paste0("U#", sub(paste0("^", strip), "", vector[i]))
    }
  }
  
  return(result_vector)
}
