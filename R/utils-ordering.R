#' sort_idents
#'
#' Sort cell identifiers (cluster IDs, cell types) with optional reference order.
#' 
#' This function sorts identifiers in two modes:
#' 
#' 1. **Simple mode** (primary_order = NULL): Sorts by numeric or alphabetic suffix after
#'    the last underscore. Purely numeric suffixes are sorted numerically; alphabetic
#'    suffixes are sorted alphabetically.
#'    
#' 2. **Reference mode** (primary_order specified): Sorts according to a reference order
#'    with special handling for identifiers containing underscores and numeric suffixes.
#'    Produces intuitive orderings for complex cell type naming schemes (e.g., "IT_A_1",
#'    "IT_A_10", "L5PT_2").
#' 
#' @param idents Character or numeric vector of identifiers to sort (e.g., cluster IDs like 
#'   c("1", "10", "2"), or cell types like c("IT_A", "L5PT", "IT_B_1")).
#' @param primary_order Character or numeric vector specifying the desired order of primary
#'   identifier components. If NULL (default), uses simple suffix-based sorting. If provided,
#'   identifiers are first matched against this order, then secondary sorting is applied
#'   based on suffixes.
#'   
#' @details
#' **Simple mode** (primary_order = NULL):
#' - Extracts suffix after last underscore
#' - Sorts numerically if all suffixes are numeric, alphabetically otherwise
#' - Returns a factor with reversed levels for compatibility with plotting functions
#' 
#' **Reference mode** (primary_order specified):
#' 
#' The sorting algorithm works in three steps:
#' 
#' 1. **Primary matching**: Each identifier is matched against the `primary_order` vector.
#'    To handle cases where shorter patterns might match first (e.g., "1" matching before "10"),
#'    the function sorts `primary_order` by decreasing string length before pattern matching.
#'    
#' 2. **Suffix extraction**: For identifiers with underscores (e.g., "IT_A_1"), the suffix
#'    after the last underscore is extracted. Numeric suffixes are converted to numbers for
#'    proper numerical sorting. Non-numeric suffixes are prefixed with "Z" to sort after
#'    numeric suffixes.
#'    
#' 3. **Combined sorting**: Identifiers are sorted first by their position in `primary_order`,
#'    then by numeric suffix value, then by alphabetic suffix value.
#' 
#' @return Character vector of sorted identifiers (reference mode) or factor with reversed
#'   levels (simple mode).
#' 
#' @export
#' @family ordering
#' 
#' @examples
#' \dontrun{
#'  # Simple mode - sort by suffix only
#'  clusters <- c("cluster_10", "cluster_1", "cluster_2")
#'  sort_idents(clusters)
#'  # Returns factor: "cluster_1", "cluster_2", "cluster_10"
#'  
#'  # Reference mode - sort purely numeric cluster IDs
#'  clusters <- c("10", "1", "2", "20")
#'  sort_idents(clusters, primary_order = c(1, 2, 10, 20))
#'  # Returns: "1", "2", "10", "20"
#'  
#'  # Reference mode - sort cell types with primary order preference
#'  celltypes <- c("L5PT_2", "IT_A_10", "IT_A_1", "L5PT_1")
#'  sort_idents(celltypes, primary_order = c("IT_A", "L5PT"))
#'  # Returns: "IT_A_1", "IT_A_10", "L5PT_1", "L5PT_2"
#'  
#'  # Reference mode - sort mixed numeric and labeled identifiers
#'  mixed <- c("10", "1", "T2", "T1")
#'  sort_idents(mixed, primary_order = c(1, 10, "T1", "T2"))
#'  # Returns: "1", "10", "T1", "T2"
#' }
sort_idents <- function(idents, primary_order = NULL) {
  
  # Simple mode: sort by suffix only
  if (is.null(primary_order)) {
    suffixes <- sub(".*_", "", idents)
    
    # Determine if suffixes are numeric or alphabetic
    if (all(grepl("^[0-9]+$", suffixes))) {
      sorted_vec <- idents[order(as.numeric(suffixes))]
    } else {
      sorted_vec <- idents[order(suffixes)]
    }
    
    return(factor(rev(sorted_vec), levels = rev(sorted_vec)))
  }
  
  # Reference mode: sort by primary_order with suffix handling
  idents <- as.character(idents)
  primary_order <- as.character(primary_order)
  
  # Sort primary_order by decreasing length to match longer patterns first
  primary_order_sorted <- primary_order[order(nchar(primary_order), decreasing = TRUE)]
  
  # Extract the primary identifier component by matching against primary_order
  primary <- sapply(idents, function(x) str_extract(x, paste(primary_order_sorted, collapse = "|")))
  
  # Extract suffix after the last underscore (if present)
  suffix <- sapply(idents, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
  
  # Convert numeric suffixes to numbers for proper numerical sorting
  suffix_numeric <- suppressWarnings(as.numeric(suffix))
  
  # For non-numeric suffixes, prefix with "Z" to sort them after numeric suffixes
  suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])
  suffix_numeric[is.na(suffix_numeric)] <- Inf
  
  # Create data frame for sorting
  df <- data.frame(
    ident = idents,
    primary = primary,
    suffix = suffix,
    suffix_numeric = suffix_numeric,
    stringsAsFactors = FALSE
  )
  
  # Sort by: (1) position in primary_order, (2) numeric suffix, (3) alphabetic suffix
  df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
  
  return(unlist(df$ident))
}


#' sort_by_reference
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
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
