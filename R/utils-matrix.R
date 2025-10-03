#' align_matrix
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
#' @param mat (auto) parameter
#' @param row_names (auto) parameter
#' @param col_names (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
align_matrix <- function(mat, row_names, col_names) {
  aligned_mat <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                        dimnames = list(row_names, col_names))
  mat_rownames <- rownames(mat)
  mat_colnames <- colnames(mat)
  for (r in mat_rownames) {
    for (c in mat_colnames) {
      aligned_mat[r, c] <- mat[r, c]
    }
  }
  return(aligned_mat)
}
