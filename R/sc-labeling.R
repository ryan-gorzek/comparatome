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
#' Propagate labels to unlabeled cells based on nearest neighbor voting in reduced dimensional space.
#' 
#' For cells marked as "None" in the specified identity column, this function finds their 
#' k nearest neighbors in a chosen embedding space (PCA, UMAP, etc.) and assigns labels 
#' if a sufficient fraction of neighbors agree on the label. This approach is useful for
#' refining cell type assignments after initial clustering or for labeling ambiguous cells
#' based on their similarity to confidently labeled cells.
#'
#' @param obj Seurat object with the specified reduction and existing labels in metadata.
#' @param ident Character string specifying the name of the metadata column containing labels.
#'   Cells marked as "None" or NA are considered unlabeled and candidates for label propagation.
#'   The original column is not modified.
#' @param output_col Character string specifying the name of the output metadata column.
#'   Default is `[ident]_nn`. Choose descriptive names when testing multiple parameter
#'   combinations (e.g., "subclass_nn_pca30", "subclass_nn_umap").
#' @param reduction Character string specifying which dimensional reduction to use for
#'   nearest neighbor search. Default is "pca". Common options include "pca", "umap",
#'   "tsne", or any other reduction stored in the Seurat object.
#' @param dims Integer vector specifying which dimensions to use from the reduction.
#'   Default is 1:30 (first 30 dimensions). For PCA, using more dimensions captures
#'   more variation but may include noise. For UMAP/tSNE, typically all dimensions
#'   (e.g., 1:2) are used.
#' @param fraction Numeric value between 0 and 1 specifying the minimum fraction of 
#'   neighbors required to agree for label assignment. Default is 0.6 (60% agreement).
#'   Higher values increase stringency and reduce false assignments.
#' @param n.neighbors Integer specifying the number of nearest neighbors to consider 
#'   for voting. Default is 20. Larger values smooth over local noise but may blur
#'   boundaries between cell types.
#'
#' @return Seurat object with new metadata column named according to `output_col`.
#'   The original `ident` column is not modified.
#'   Output column values:
#'   \itemize{
#'     \item Originally labeled cells: NA (no change needed)
#'     \item Unlabeled cells meeting threshold: assigned label
#'     \item Unlabeled cells not meeting threshold: "None"
#'   }
#'
#' @details
#' **Label propagation workflow:**
#' 
#' 1. Extract coordinates from specified reduction (e.g., PC1-PC30)
#' 2. Identify labeled cells (ident != "None" and !is.na(ident)) vs unlabeled cells
#' 3. For each unlabeled cell, find k nearest neighbors using FNN::get.knnx
#' 4. Calculate the fraction of neighbors with each label
#' 5. Assign the most common label if it exceeds the threshold fraction
#' 6. Otherwise, keep the cell labeled as "None"
#' 
#' **Choosing the reduction and dimensions:**
#' 
#' - **PCA** (default, dims = 1:30): Recommended for most cases. PCA captures global
#'   structure and is less sensitive to local noise. Using 20-50 dimensions balances
#'   biological signal with noise. PCA-based labeling works well even when UMAP/tSNE
#'   show complex structure.
#'   
#' - **UMAP** (dims = 1:2): Uses the 2D UMAP embedding. Emphasizes local neighborhood
#'   structure and can be effective when cell types form distinct visual clusters.
#'   However, UMAP can create artificial separation, so use with caution.
#'   
#' - **Custom reductions**: Any reduction in the Seurat object can be used. For
#'   integrated analyses, "integrated" or "harmony" reductions are appropriate.
#' 
#' **Parameter tuning:**
#' 
#' - Increase `fraction` (e.g., 0.7-0.8) for more conservative labeling
#' - Increase `n.neighbors` (e.g., 30-50) for smoother boundaries
#' - Decrease both for more aggressive labeling of ambiguous cells
#' 
#' **Limitations:**
#' 
#' - Assumes labeled cells are representative of all cell states
#' - Cannot identify novel cell types (will assign to nearest existing type)
#' - Performance depends on quality of initial labels and embedding
#' 
#' @export
#' @family labeling
#'
#' @examples
#' \dontrun{
#'   # Default: Use PCA (first 30 dimensions) for label propagation
#'   obj <- LabelByNearestNeighbors(
#'     obj, 
#'     ident = "subclass",
#'     output_col = "subclass_nn",
#'     fraction = 0.7,
#'     n.neighbors = 30
#'   )
#'   
#'   # Original column unchanged
#'   any(is.na(obj$subclass))  # Still TRUE if there were NAs
#'   
#'   # Test with UMAP embedding
#'   obj <- LabelByNearestNeighbors(
#'     obj,
#'     ident = "subclass",
#'     output_col = "subclass_nn_umap",
#'     reduction = "umap",
#'     dims = 1:2,
#'     fraction = 0.6,
#'     n.neighbors = 20
#'   )
#'   
#'   # Compare results
#'   table(obj$subclass_nn, useNA = "ifany")
#'   table(obj$subclass_nn_umap, useNA = "ifany")
#' }
LabelByNearestNeighbors <- function(obj, ident, output_col = NULL, reduction = "pca", 
                                    dims = 1:30, fraction = 0.6, n.neighbors = 20) {
  
  library(Seurat)
  library(FNN)
  library(dplyr)
  
  # Set default output column name if not provided
  if (is.null(output_col)) {
    output_col <- paste0(ident, "_nn")
  }
  
  # Work with a copy of the labels to avoid modifying original column
  labels <- obj[[ident]][, 1]
  
  # Treat NA and "None" as unlabeled
  labels[is.na(labels)] <- "None"
  
  # Validate that the reduction exists
  if (!reduction %in% names(obj@reductions)) {
    stop(paste0("Reduction '", reduction, "' not found in Seurat object. ",
                "Available reductions: ", paste(names(obj@reductions), collapse = ", ")))
  }
  
  # Extract embeddings from specified reduction
  embeddings <- Embeddings(obj, reduction)
  
  # Validate dimensions
  max_dim <- ncol(embeddings)
  if (max(dims) > max_dim) {
    stop(paste0("Requested dimensions exceed available dimensions. ",
                "Reduction '", reduction, "' has ", max_dim, " dimensions, ",
                "but dims requests up to ", max(dims)))
  }
  
  # Subset to requested dimensions
  embeddings <- embeddings[, dims, drop = FALSE]
  
  # Identify labeled and unlabeled cells using our working copy
  labeled_cells <- rownames(obj@meta.data)[labels != "None"]
  unlabeled_cells <- rownames(obj@meta.data)[labels == "None"]
  all_cells <- c(labeled_cells, unlabeled_cells)
  
  # Find k nearest neighbors for each unlabeled cell
  nn <- get.knnx(
    data = embeddings[all_cells, , drop = FALSE],      # All cells as reference
    query = embeddings[unlabeled_cells, , drop = FALSE], # Unlabeled cells as query
    k = n.neighbors
  )
  
  # For each unlabeled cell, collect labels from neighbors
  neighbor_idents <- vector("list", length = nrow(nn$nn.index))
  
  for (i in 1:nrow(nn$nn.index)) {
    neighbor_indices <- nn$nn.index[i, ]
    neighbor_cells <- all_cells[neighbor_indices]
    # Use working copy of labels
    neighbor_idents[[i]] <- labels[match(neighbor_cells, rownames(obj@meta.data))]
  }
  
  # Calculate fraction of neighbors with each label
  all_idents <- unique(labels[match(all_cells, rownames(obj@meta.data))])
  neighbor_subclass_fractions <- lapply(neighbor_idents, function(idents) {
    tbl <- table(factor(idents, levels = all_idents))
    prop.table(tbl)
  })
  
  # Combine into dataframe
  neighbor_fraction_df <- do.call(rbind, neighbor_subclass_fractions)
  rownames(neighbor_fraction_df) <- unlabeled_cells
  
  # Assign labels based on threshold
  assigned_idents <- sapply(1:nrow(neighbor_fraction_df), function(i) {
    fractions <- neighbor_fraction_df[i, ]
    fractions_no_none <- fractions[names(fractions) != "None"]
    
    # If no labeled neighbors, keep as "None"
    if (length(fractions_no_none) == 0 || all(fractions_no_none == 0)) {
      return("None")
    }
    
    # Find most common label
    top_subclass <- names(fractions_no_none)[which.max(fractions_no_none)]
    
    # Assign if threshold is met
    if (max(fractions_no_none) >= fraction) {
      return(top_subclass)
    } else {
      return("None")
    }
  })
  
  # Add results to metadata
  # Originally labeled cells get NA, unlabeled cells get assigned labels or "None"
  obj[[output_col]] <- "None"
  obj[[output_col]][rownames(neighbor_fraction_df), 1] <- assigned_idents
  obj[[output_col]][labeled_cells, 1] <- NA
  
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
