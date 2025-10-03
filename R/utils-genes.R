#' top_genes_two
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
#' @param df (auto) parameter
#' @param sort_column_asc (auto) parameter
#' @param sort_column_desc (auto) parameter
#' @param idx (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
top_genes_two <- function(df, sort_column_asc, sort_column_desc, idx) {
  
  # Sort the dataframe by the specified columns: one ascending and one descending
  sorted_df <- df[order(df[[sort_column_asc]], -df[[sort_column_desc]]), ]
  
  # Get the top n rownames
  top_n <- sorted_df$gene[idx]
  
  return(top_n[!is.na(top_n)])
}


#' top_genes_asc
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
#' @param df (auto) parameter
#' @param sort_column_asc (auto) parameter
#' @param idx (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
top_genes_asc <- function(df, sort_column_asc, idx) {
  
  # Sort the dataframe by the specified columns: one ascending and one descending
  sorted_df <- df[order(df[[sort_column_asc]]), ]
  
  # Get the top n rownames
  top_n <- sorted_df$gene[idx]
  
  return(top_n[!is.na(top_n)])
}


#' top_genes_desc
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
#' @param df (auto) parameter
#' @param sort_column_desc (auto) parameter
#' @param idx (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
top_genes_desc <- function(df, sort_column_desc, idx) {
  
  # Sort the dataframe by the specified columns: one ascending and one descending
  sorted_df <- df[order(-df[[sort_column_desc]]), ]
  
  # Get the top n rownames
  top_n <- sorted_df$gene[idx]
  
  return(top_n[!is.na(top_n)])
}
