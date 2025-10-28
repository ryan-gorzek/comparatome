#' sort_idents
#'
#' Sort cell identifiers (cluster IDs, cell types) according to a specified reference order, 
#' handling mixed numeric and alphanumeric patterns.
#' 
#' This function sorts identifiers by matching them against a reference order (primary_order),
#' with special handling for identifiers containing underscores and numeric suffixes. It is
#' designed to produce intuitive orderings for complex cell type naming schemes commonly used
#' in single-cell analyses (e.g., "IT_A_1", "IT_A_10", "L5PT_2").
#' 
#' @param idents Character or numeric vector of identifiers to sort (e.g., cluster IDs like 
#'   c("1", "10", "2"), or cell types like c("IT_A", "L5PT", "IT_B_1")).
#' @param primary_order Character or numeric vector specifying the desired order of primary
#'   identifier components. Identifiers are first matched against this order, then secondary
#'   sorting is applied based on suffixes.
#'   
#' @details
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
#' **Examples of sorting behavior**:
#' - Purely numeric: c("2", "1", "10") becomes c("1", "2", "10")
#' - Mixed with prefixes: c("IT_2", "IT_1", "IT_10") becomes c("IT_1", "IT_2", "IT_10")
#' - Multiple prefixes: c("L5PT_1", "IT_A_1", "IT_B_1") sorted according to primary_order
#' 
#' @return Character vector of sorted identifiers in the order determined by the sorting algorithm.
#' 
#' @export
#' @family ordering
#' 
#' @examples
#' \dontrun{
#'  # Sort purely numeric cluster IDs
#'  clusters <- c("10", "1", "2", "20")
#'  sort_idents(clusters, primary_order = c(1, 2, 10, 20))
#'  # Returns: "1", "2", "10", "20"
#'  
#'  # Sort cell types with primary order preference
#'  celltypes <- c("L5PT_2", "IT_A_10", "IT_A_1", "L5PT_1")
#'  sort_idents(celltypes, primary_order = c("IT_A", "L5PT"))
#'  # Returns: "IT_A_1", "IT_A_10", "L5PT_1", "L5PT_2"
#'  
#'  # Sort mixed numeric and labeled identifiers
#'  mixed <- c("10", "1", "T2", "T1")
#'  sort_idents(mixed, primary_order = c(1, 10, "T1", "T2"))
#'  # Returns: "1", "10", "T1", "T2"
#' }
sort_idents <- function(idents, primary_order) {
  
  # Convert everything to character for consistent handling
  idents <- as.character(idents)
  primary_order <- as.character(primary_order)
  
  # Sort primary_order by decreasing length to match longer patterns first
  # This prevents "1" from matching before "10" in regex patterns
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


#' #' sort_idents
#' #'
#' #' Auto-generated roxygen skeleton for comparatome.
#' #' Part of the helpers family.
#' #' @param vec (auto) parameter
#' #' @return (auto) value; see function body.
#' #' @keywords internal
#' #' @family helpers
#' #' @examples
#' #' \dontrun{
#' #'  # Example usage will be added
#' #' }
#' sort_idents <- function(vec) {
#'   suffixes <<- sub(".*_", "", vec)
#'   
#'   # Determine if suffixes are numeric or alphabetic
#'   if (all(grepl("^[0-9]+$", suffixes))) {
#'     sorted_vec <<- vec[order(as.numeric(suffixes))]
#'   } else {
#'     sorted_vec <<- vec[order(suffixes)]
#'   }
#'   
#'   return(factor(rev(sorted_vec), levels = rev(sorted_vec)))
#' }


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
