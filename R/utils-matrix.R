#' add_zeros_to_table
#'
#' Expand a matrix or table by adding zero-filled rows and columns for missing labels.
#' 
#' This utility function ensures a matrix contains all specified row and column names,
#' adding zero-filled rows/columns for any missing names. This is particularly useful
#' for confusion matrices where not all classes may be represented in predictions,
#' but a complete matrix is needed for visualization or downstream analysis.
#' 
#' @param tbl Matrix or table to expand
#' @param new_row_names Character vector of row names that should exist in the table.
#'   Any names not already present will be added as zero-filled rows.
#' @param new_col_names Character vector of column names that should exist in the table.
#'   Any names not already present will be added as zero-filled columns.
#'   
#' @details
#' The function operates in two phases:
#' 
#' 1. **Row expansion**: For each name in `new_row_names` not present in the table,
#'    adds a new row filled with zeros. Row names are preserved and the new row name
#'    is appended.
#'    
#' 2. **Column expansion**: For each name in `new_col_names` not present in the table,
#'    adds a new column filled with zeros. Column names are preserved and the new
#'    column name is appended.
#' 
#' The order of rows and columns in the output follows the original order, with new
#' rows/columns appended at the end. Use standard matrix subsetting to reorder if needed.
#' 
#' @return Matrix with the same structure as `tbl` but expanded to include all specified
#'   row and column names, with zeros filling any newly added cells.
#'   
#' @keywords internal
#' @family matrix
#' 
#' @examples
#' \dontrun{
#'  # Create a confusion matrix missing some classes
#'  conf_matrix <- matrix(c(10, 2, 3, 15), nrow = 2,
#'                        dimnames = list(c("A", "B"), c("X", "Y")))
#'  
#'  # Expand to include all expected classes
#'  complete_matrix <- add_zeros_to_table(
#'    conf_matrix,
#'    new_row_names = c("A", "B", "C"),
#'    new_col_names = c("X", "Y", "Z")
#'  )
#'  # Result includes zero-filled row for "C" and column for "Z"
#'  
#'  # Common use case: ensure mapping heatmap shows all reference types
#'  mapping_table <- table(query_labels, predicted_labels)
#'  complete_table <- add_zeros_to_table(
#'    mapping_table,
#'    new_row_names = all_query_types,
#'    new_col_names = all_reference_types
#'  )
#' }
add_zeros_to_table <- function(tbl, new_row_names, new_col_names) {
  
  # Add new rows of zeros for any missing row names
  for (row_name in new_row_names) {
    if (!any(row_name %in% rownames(tbl))) {
      # Create zero-filled row with same number of columns as table
      tbl <- rbind(tbl, setNames(t(rep(0, ncol(tbl))), row_name))
      # Preserve original row names and append new name
      orig_names <- rownames(tbl)[rownames(tbl) != ""]
      rownames(tbl) <- c(orig_names, row_name)
    }
  }
  
  # Add new columns of zeros for any missing column names
  for (col_name in new_col_names) {
    if (!any(col_name %in% colnames(tbl))) {
      # Create zero-filled column with same number of rows as table
      tbl <- cbind(tbl, setNames(rep(0, nrow(tbl)), col_name))
      # Preserve original column names and append new name
      orig_names <- colnames(tbl)[colnames(tbl) != ""]
      colnames(tbl) <- c(orig_names, col_name)
    }
  }
  
  return(tbl)
}


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
