#' LabelCells
#'
#' Assign final subclass and type labels to cells based on their clustering at specified resolutions.
#' For each subclass, cells are labeled with that subclass name, and types are assigned by ordering
#' clusters by size (largest cluster becomes [subclass]_1, second largest becomes [subclass]_2, etc.).
#'
#' @param obj Seurat object with metadata columns formatted as "subclass.[resolution]" and "SCT_snn_res.[resolution]"
#' @param subclass_resolution Named list where names are subclass labels and values are clustering resolutions 
#'   (e.g., list(IT_A = 0.5, L5PT = 1.5)). The resolution indicates which clustering result to use for type assignment.
#'
#' @return Seurat object with three new metadata columns:
#' \itemize{
#'   \item subclass: The subclass label (e.g., "IT_A", "Pvalb")
#'   \item type: Subclass + cluster rank (e.g., "IT_A_1", "Pvalb_2")
#'   \item subclass.type: Copy of subclass label
#' }
#'
#' @details
#' The function operates by:
#' \enumerate{
#'   \item For each subclass-resolution pair, identifying cells with that subclass assignment
#'   \item Extracting cluster identities from the specified resolution
#'   \item Ranking clusters by cell count (descending)
#'   \item Assigning type labels as [subclass]_[rank]
#' }
#'
#' If a subclass contains only one cluster, type label equals subclass label.
#' Issues a warning if any NAs are present in the final type column.
#'
#' @export
#' @family labeling
#'
#' @examples
#' \dontrun{
#'   # After running SubclassByIdent to assign subclasses at multiple resolutions
#'   subclass.resolutions <- rev(list(
#'     IT_A = 0.5,
#'     IT_B = 0.5,
#'     L5PT = 1.5,
#'     L6CT = 1.5
#'   ))
#'   obj <- LabelCells(obj, subclass.resolutions)
#'   
#'   # Check results
#'   table(obj$subclass, obj$type)
#' }
LabelCells <- function(obj, subclass_resolution) {

  # Iterate over each subclass and resolution
  for (subclass in names(subclass_resolution)) {
    resolution <- subclass_resolution[[subclass]]
    
    # Extract the subclass and cluster columns based on the resolution
    subclass_col <- paste0("subclass.", resolution)
    cluster_col <- paste0("SCT_snn_res.", resolution)
    
    # Get the cells belonging to the current subclass
    cells <- obj@meta.data[[subclass_col]] == subclass
    
    # Assign the subclass label to the 'subclass' column
    obj@meta.data$subclass[cells] <- subclass
    
    # Extract the cluster labels for these cells
    clusters <- obj@meta.data[[cluster_col]][cells]
    
    # Get the cluster sizes in decreasing order
    cluster_sizes <- sort(table(clusters), decreasing = TRUE)
    
    # Create a mapping from original cluster labels to new type labels
    cluster_to_type <- setNames(paste0(subclass, "_", seq_along(cluster_sizes)), names(cluster_sizes))
    
    # Assign the new type labels based on the cluster sizes
    if (sum(cluster_sizes > 0) > 1) {
      obj@meta.data$type[cells] <- cluster_to_type[as.character(clusters)]
    } else { obj@meta.data$type[cells] <- subclass }
    obj@meta.data$subclass.type[cells] <- subclass
  }
  
  if (any(is.na(obj@meta.data$type))) { warning("NAs present in type column...") }
  return(obj)

}


#' LabelByNearestNeighbors
#'
#' Propagate labels to unlabeled cells based on UMAP nearest neighbor voting.
#' For cells marked as "None" in the specified identity column, this function finds their 
#' k nearest neighbors in UMAP space and assigns labels if a sufficient fraction of 
#' neighbors agree on the label.
#'
#' @param obj Seurat object with UMAP reduction and existing labels
#' @param ident Character, name of metadata column containing labels (with some cells marked "None")
#' @param fraction Numeric, minimum fraction of neighbors required to agree for label assignment (default: 0.6)
#' @param n.neighbors Integer, number of nearest neighbors to consider for voting (default: 20)
#'
#' @return Seurat object with new metadata column named [ident]_nn containing propagated labels.
#'   Cells that were originally labeled remain NA in the new column.
#'   Cells that don't meet the threshold remain "None".
#'
#' @details
#' Label propagation workflow:
#' \enumerate{
#'   \item Extract UMAP coordinates for all cells
#'   \item Identify labeled vs unlabeled cells (unlabeled = "None")
#'   \item For each unlabeled cell, find k nearest neighbors using FNN::get.knnx
#'   \item Calculate the fraction of neighbors with each label
#'   \item If the top label exceeds the threshold, assign it; otherwise keep "None"
#' }
#'
#' Requires the FNN package for efficient nearest neighbor search.
#' Originally labeled cells receive NA in the output column.
#'
#' @export
#' @family labeling
#'
#' @examples
#' \dontrun{
#'   # After partial labeling, propagate to ambiguous cells
#'   obj <- LabelByNearestNeighbors(
#'     obj, 
#'     ident = "subclass",
#'     fraction = 0.7,  # Require 70% agreement
#'     n.neighbors = 30
#'   )
#'   
#'   # Check how many cells were labeled
#'   table(obj$subclass_nn)
#' }
LabelByNearestNeighbors <- function(obj, ident, fraction = 0.6, n.neighbors = 20) {
  
  library(Seurat)
  library(FNN)  # For fast nearest neighbor search
  library(dplyr)
  
  obj[[ident]][is.na(obj[[ident]])] <- "None"
  
  # 1. Get UMAP embeddings
  umap_embeddings <- Embeddings(obj, "umap")
  
  # Extract metadata
  meta <- obj@meta.data
  
  # Get cell names directly
  labeled_cells <- rownames(meta)[meta[[ident]] != "None"]
  unlabeled_cells <- rownames(meta)[meta[[ident]] == "None"]
  all_cells <- c(labeled_cells, unlabeled_cells)
  
  # 3. Find k nearest neighbors for each unlabeled cell
  k <- n.neighbors  # Number of neighbors to consider
  # Nearest neighbors: query unlabeled cells against labeled cells
  nn <- get.knnx(
    data = umap_embeddings[all_cells, ],    # Reference (labeled cells)
    query = umap_embeddings[unlabeled_cells, ], # Query (unlabeled cells)
    k = k
  )
  
  # 4. For each unlabeled cell, check the subclass of neighbors
  neighbor_idents <- vector("list", length = nrow(nn$nn.index))
  
  for (i in 1:nrow(nn$nn.index)) {
    neighbor_indices <- nn$nn.index[i, ]
    neighbor_cells <- all_cells[neighbor_indices]
    
    # Save subclass labels
    neighbor_idents[[i]] <- obj[[ident]][neighbor_cells, 1]
  }
  
  # 5. Summarize neighbor subclass proportions per cell
  all_idents <- unique(obj[[ident]][all_cells, 1])
  neighbor_subclass_fractions <- lapply(neighbor_idents, function(idents) {
    tbl <- table(factor(idents, levels = all_idents))
    prop.table(tbl)
  })
  
  # Combine into a dataframe
  neighbor_fraction_df <- do.call(rbind, neighbor_subclass_fractions)
  rownames(neighbor_fraction_df) <- unlabeled_cells
  
  # Threshold for assignment
  threshold <- fraction
  
  assigned_idents <- sapply(1:nrow(neighbor_fraction_df), function(i) {
    fractions <- neighbor_fraction_df[i, ]
    fractions_no_none <- fractions[names(fractions) != "None"]
    
    if (length(fractions_no_none) == 0 || all(fractions_no_none == 0)) {
      return("None")
    }
    
    top_subclass <- names(fractions_no_none)[which.max(fractions_no_none)]
    
    if (max(fractions_no_none) >= threshold) {
      return(top_subclass)
    } else {
      return("None")
    }
  })
  
  # Add back to Seurat object's metadata
  obj[[paste0(ident, "_nn")]] <- "None"
  obj[[paste0(ident, "_nn")]][rownames(neighbor_fraction_df), 1] <- assigned_idents
  obj[[paste0(ident, "_nn")]][labeled_cells, 1] <- NA
  
  return(obj)
}


#' IdentBySample
#'
#' Plot relative proportions of each identity (cluster/subclass/type) across samples.
#' Creates a bar plot showing median proportions with individual sample values as scatter points.
#'
#' @param obj Seurat object with 'sample' metadata and active identity set
#' @param y_limits Numeric vector of length 2 specifying y-axis limits (default: c(0, 0.50))
#'
#' @return NULL (plots are printed directly). Generates a ggplot2 bar plot with:
#' \itemize{
#'   \item Bars showing median proportion of each identity across samples
#'   \item Points showing individual sample values
#'   \item Percentage labels above bars
#' }
#'
#' @details
#' For each identity level and sample, calculates the proportion of cells.
#' Then computes median proportion across samples for plotting.
#' Drops cells with NA in the active identity and warns if any are found.
#' Identity order is taken from the Seurat object's levels.
#'
#' @export
#' @family labeling
#'
#' @examples
#' \dontrun{
#'   # Plot subclass proportions
#'   Idents(obj) <- "subclass"
#'   levels(obj) <- c("IT_A", "IT_B", "L5PT", "L6CT")
#'   IdentBySample(obj, y_limits = c(0, 0.60))
#' }
IdentBySample <- function(obj, y_limits = c(0, 0.50)) {
  
  # Assuming your dataframe is named df with columns 'subclass' and 'sample'
  df <- obj[[]]
  df$active.ident <- obj@active.ident
  df$smpl <- df$sample
  specified_order <- levels(obj)

  # Check for NAs and warn the user
  na_samples <- df %>% filter(is.na(active.ident)) %>% dplyr::count(smpl)
  if(nrow(na_samples) > 0) {
    warning("The following samples contained NAs and were dropped: ",
            paste(na_samples$sample, " (", na_samples$n, ")", collapse = "\n"))
  }
  
  # Drop NAs
  df <- df %>% drop_na(active.ident)
  
  # Calculate the relative proportions of each active.ident by sample
  relative_proportions <- df %>%
    group_by(sample, active.ident) %>%
    dplyr::summarize(count = n(), .groups = 'drop') %>%
    ungroup() %>%
    group_by(sample) %>%
    mutate(Proportion = count / sum(count))
  
  # Calculate the median proportions across samples
  median_proportions <- relative_proportions %>%
    group_by(active.ident) %>%
    dplyr::summarize(MedianProportion = median(Proportion))
  
  # Rename the columns for clarity
  colnames(median_proportions) <- c("active.ident", "MedianProportion")
  
  # Convert active.ident to a factor and specify the order
  median_proportions$active.ident <- factor(median_proportions$active.ident, levels = specified_order)
  relative_proportions$active.ident <- factor(relative_proportions$active.ident, levels = specified_order)
  
  # Calculate the maximum proportion for each active.ident
  max_proportions <- relative_proportions %>%
    group_by(active.ident) %>%
    summarize(MaxProportion = max(Proportion))
  
  # Merge the max proportions with the median proportions
  plot_data <- merge(median_proportions, max_proportions, by = "active.ident")
  
  # Create the barplot
  p <- ggplot(plot_data, aes(x = active.ident, y = MedianProportion, fill = active.ident)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = scales::percent(MedianProportion, accuracy = 0.1), 
                  y = MaxProportion + 0.01), # Place text labels slightly above the highest scatter point
              vjust = -0.5) + # Adjust the position of the text
    geom_point(data = relative_proportions, aes(x = active.ident, y = Proportion), 
               color = "black") + # Scatter the different samples
    theme_minimal() +
    theme(legend.position = "none", 
          aspect.ratio = 1) + # Remove the legend and make the plot square
    xlab("") +
    ylab("Relative Proportion") +
    ggtitle("") +
    scale_y_continuous(labels = scales::percent, limits = y_limits) +
    scale_x_discrete(limits = specified_order) + 
    theme(axis.text.x = element_text(size=12), 
          axis.text.y = element_text(size=12), 
          axis.title.x = element_text(size=12), 
          axis.title.y = element_text(size=12))
  
  print(p)
  
}


#' SubclassByIdent
#'
#' Assign subclass labels to cells based on their cluster membership at specified resolutions.
#' Creates metadata columns formatted as "subclass.[resolution]" containing subclass assignments.
#'
#' @param obj Seurat object with clustering results stored as "SCT_snn_res.[resolution]"
#' @param subclass.idx Nested list structure:
#'   List level 1: resolution names (e.g., "SCT_snn_res.0.2")
#'   List level 2: subclass names (e.g., "IT_A", "Pvalb")
#'   Values: character vectors of cluster IDs belonging to that subclass
#'
#' @return Seurat object with new metadata columns named "subclass.[resolution]".
#'   Cells in specified clusters receive the subclass label; others remain NA.
#'
#' @details
#' Example structure for subclass.idx:
#' \preformatted{
#' subclass.idx <- list(
#'   SCT_snn_res.0.2 = list(
#'     IT_A = c("1", "2"),
#'     L5PT = c("6")
#'   ),
#'   SCT_snn_res.0.5 = list(
#'     IT_A = c("2", "3", "4"),
#'     L5PT = c("9")
#'   )
#' )
#' }
#'
#' The function automatically creates appropriately named metadata columns based on
#' the resolution format (handles both "res.X" and "res.X.Y" patterns).
#'
#' @export
#' @family labeling
#'
#' @examples
#' \dontrun{
#'   subclass.idx <- list()
#'   subclass.idx$SCT_snn_res.0.2$IT_A <- c("1", "2")
#'   subclass.idx$SCT_snn_res.0.2$L5PT <- c("6")
#'   subclass.idx$SCT_snn_res.0.5$IT_A <- c("2", "3", "4", "6")
#'   subclass.idx$SCT_snn_res.0.5$L5PT <- c("9")
#'   
#'   obj <- SubclassByIdent(obj, subclass.idx)
#'   
#'   # Check assignments at one resolution
#'   table(obj$subclass.0.2)
#' }
SubclassByIdent <- function(obj, subclass.idx) {
  
  for (id in names(subclass.idx)) {
    
    Idents(obj) <- id
    
    id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
    if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
    else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
    else { sbcl.col <- paste0("subclass.", id.num[1]) }
    
    if ((sbcl.col %in% colnames(obj[[]]) == FALSE)) { obj[[sbcl.col]] <- NA }
    
    for (sbcl in names(subclass.idx[[id]])) {
      
      if (!is.null(subclass.idx[[id]][[sbcl]])) {
        
        cell.names <- WhichCells(obj, ident = subclass.idx[[id]][[sbcl]])
        obj[[sbcl.col]][cell.names,] <- sbcl
        
      }
    }
  }
  
  return(obj)
  
}
