#' URLencodeNCBI
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
#' @param x (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
URLencodeNCBI <- function(x) {
  b <- charToRaw(enc2utf8(x))
  out <- vapply(b, function(xx) {
    if ((xx >= 0x30 && xx <= 0x39) || (xx >= 0x41 && xx <= 0x5A) || (xx >= 0x61 && xx <= 0x7A)) {
      rawToChar(as.raw(xx))
    } else {
      sprintf("%%%02X", as.integer(xx))
    }
  }, character(1))
  paste0(out, collapse = "")
}


#' make_folder
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
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
#' Part of the helpers family.
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
