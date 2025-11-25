#' top_genes_desc
#'
#' Extract top genes from a data frame sorted by a column in descending order.
#' Commonly used to select highly ranked marker genes or top differentially expressed genes.
#'
#' @param df Data frame containing gene information with a 'gene' column and sorting metric.
#'   Typically output from FindMarkers or similar differential expression functions.
#' @param sort_column_desc String specifying the column name to sort by in descending order.
#'   Common choices include "avg_log2FC" for fold change or "pct.diff" for expression difference.
#' @param idx Numeric vector of indices to extract after sorting (e.g., 1:20 for top 20 genes).
#'
#' @return Character vector of gene names corresponding to the specified indices after sorting,
#'   with NA values removed.
#'
#' @details
#' The function:
#' 1. Sorts the input data frame by the specified column in descending order
#' 2. Extracts the 'gene' column values at the requested indices
#' 3. Removes any NA values from the result
#'
#' This is particularly useful when creating marker gene plots where you want to visualize
#' the top N genes ranked by fold change or other metrics.
#'
#' @keywords internal
#' @family helpers
#'
#' @examples
#' \dontrun{
#'  # Find markers for a cluster
#'  markers <- FindMarkers(obj, ident.1 = "IT_A", ident.2 = "IT_B")
#'  markers$gene <- rownames(markers)
#'  markers$pct.diff <- markers$pct.1 - markers$pct.2
#'  
#'  # Get top 20 genes by fold change
#'  top_fc <- top_genes_desc(markers, "avg_log2FC", 1:20)
#'  
#'  # Get top 10 genes by expression difference
#'  top_pct <- top_genes_desc(markers, "pct.diff", 1:10)
#' }
top_genes_desc <- function(df, sort_column_desc, idx) {
  
  # Sort the dataframe by the specified column in descending order
  sorted_df <- df[order(-df[[sort_column_desc]]), ]
  
  # Get the top n gene names at the specified indices
  top_n <- sorted_df$gene[idx]
  
  # Remove NA values and return
  return(top_n[!is.na(top_n)])
}


#' top_genes_asc
#'
#' Extract top genes from a data frame sorted by a column in ascending order.
#' Useful for selecting genes with lowest p-values or other ascending metrics.
#'
#' @param df Data frame containing gene information with a 'gene' column and sorting metric.
#'   Typically output from FindMarkers or similar differential expression functions.
#' @param sort_column_asc String specifying the column name to sort by in ascending order.
#'   Common choices include "p_val_adj" for adjusted p-values.
#' @param idx Numeric vector of indices to extract after sorting (e.g., 1:20 for top 20 genes).
#'
#' @return Character vector of gene names corresponding to the specified indices after sorting,
#'   with NA values removed.
#'
#' @details
#' The function:
#' 1. Sorts the input data frame by the specified column in ascending order
#' 2. Extracts the 'gene' column values at the requested indices
#' 3. Removes any NA values from the result
#'
#' This is the ascending complement to top_genes_desc, useful when lower values indicate
#' higher importance (e.g., p-values, FDR).
#'
#' @keywords internal
#' @family helpers
#'
#' @examples
#' \dontrun{
#'  # Find markers for a cluster
#'  markers <- FindMarkers(obj, ident.1 = "Pvalb", ident.2 = "Sst")
#'  markers$gene <- rownames(markers)
#'  
#'  # Get top 20 genes by adjusted p-value (most significant)
#'  top_sig <- top_genes_asc(markers, "p_val_adj", 1:20)
#' }
top_genes_asc <- function(df, sort_column_asc, idx) {
  
  # Sort the dataframe by the specified column in ascending order
  sorted_df <- df[order(df[[sort_column_asc]]), ]
  
  # Get the top n gene names at the specified indices
  top_n <- sorted_df$gene[idx]
  
  # Remove NA values and return
  return(top_n[!is.na(top_n)])
}


#' top_genes_two
#'
#' Extract top genes from a data frame sorted by two columns simultaneously.
#' Primary sort is ascending, secondary sort is descending.
#'
#' @param df Data frame containing gene information with a 'gene' column and two sorting metrics.
#'   Typically output from FindMarkers or similar differential expression functions.
#' @param sort_column_asc String specifying the primary column name to sort by in ascending order.
#'   Often used for p-values or FDR to prioritize statistical significance.
#' @param sort_column_desc String specifying the secondary column name to sort by in descending order.
#'   Often used for fold change to prioritize effect size among equally significant genes.
#' @param idx Numeric vector of indices to extract after sorting (e.g., 1:20 for top 20 genes).
#'
#' @return Character vector of gene names corresponding to the specified indices after sorting,
#'   with NA values removed.
#'
#' @details
#' The function:
#' 1. Performs a two-level sort: first by ascending order of sort_column_asc,
#'    then by descending order of sort_column_desc
#' 2. Extracts the 'gene' column values at the requested indices
#' 3. Removes any NA values from the result
#'
#' This is useful for selecting genes that are both statistically significant (low p-value)
#' and have large effect sizes (high fold change), providing a balanced ranking.
#'
#' @keywords internal
#' @family helpers
#'
#' @examples
#' \dontrun{
#'  # Find markers for a cluster
#'  markers <- FindMarkers(obj, ident.1 = "L5PT", ident.2 = "L6CT")
#'  markers$gene <- rownames(markers)
#'  
#'  # Get top 20 genes prioritizing significance, then fold change
#'  top_balanced <- top_genes_two(markers, "p_val_adj", "avg_log2FC", 1:20)
#' }
top_genes_two <- function(df, sort_column_asc, sort_column_desc, idx) {
  
  # Sort the dataframe by two columns: ascending for first, descending for second
  sorted_df <- df[order(df[[sort_column_asc]], -df[[sort_column_desc]]), ]
  
  # Get the top n gene names at the specified indices
  top_n <- sorted_df$gene[idx]
  
  # Remove NA values and return
  return(top_n[!is.na(top_n)])
}
