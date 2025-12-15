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
#' @param reference_pool Character string specifying which cells to use as reference for
#'   nearest neighbor search. Options:
#'   \itemize{
#'     \item "all" (default): Use all cells (labeled + unlabeled) as reference pool.
#'       This allows unlabeled cells to find each other as neighbors. Use when cells
#'       integrate well together and you want to propagate labels across similar unlabeled cells.
#'     \item "labeled": Use only already-labeled cells as reference pool. More conservative
#'       approach that prevents poorly-integrated cells from influencing each other.
#'       Recommended for iterative integration workflows where unlabeled cells may be
#'       poorly integrated or represent novel cell types.
#'   }
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
#' 3. Determine reference pool based on `reference_pool` parameter
#' 4. For each unlabeled cell, find k nearest neighbors in reference pool
#' 5. Calculate the fraction of neighbors with each label
#' 6. Assign the most common label if it exceeds the threshold fraction
#' 7. Otherwise, keep the cell labeled as "None"
#' 
#' **Choosing reference_pool:**
#' 
#' - **"all"**: Default behavior. Uses all cells as reference pool. Appropriate when:
#'   - All cells integrate well in the embedding space
#'   - You want to propagate labels across clusters of similar unlabeled cells
#'   - Running single-pass labeling on well-integrated data
#'   
#' - **"labeled"**: Uses only labeled cells as reference. Recommended when:
#'   - Running iterative integration with multiple rounds
#'   - Unlabeled cells form separate clusters in UMAP (poor integration)
#'   - You want to prevent error propagation from ambiguous cells
#'   - Unlabeled cells may represent novel cell types or technical artifacts
#'   - Working with spatial transcriptomics where spatial cells may not integrate perfectly
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
#' - Use `reference_pool = "labeled"` for more conservative labeling in iterative workflows
#' 
#' **Limitations:**
#' 
#' - Assumes labeled cells are representative of all cell states
#' - Cannot identify novel cell types (will assign to nearest existing type or "None")
#' - Performance depends on quality of initial labels and embedding
#' - With `reference_pool = "all"`, poorly integrated cells can influence each other
#' 
#' @export
#' @family labeling
#'
#' @examples
#' \dontrun{
#'   # Default: Use all cells as reference
#'   obj <- LabelByNearestNeighbors(
#'     obj, 
#'     ident = "subclass",
#'     output_col = "subclass_nn",
#'     fraction = 0.7,
#'     n.neighbors = 30,
#'     reference_pool = "all"
#'   )
#'   
#'   # Conservative: Use only labeled cells as reference (recommended for spatial)
#'   obj <- LabelByNearestNeighbors(
#'     obj,
#'     ident = "subclass",
#'     output_col = "subclass_nn_labeled",
#'     fraction = 0.6,
#'     n.neighbors = 100,
#'     reference_pool = "labeled"
#'   )
#'   
#'   # Test with UMAP embedding
#'   obj <- LabelByNearestNeighbors(
#'     obj,
#'     ident = "subclass",
#'     output_col = "subclass_nn_umap",
#'     reduction = "umap",
#'     dims = 1:2,
#'     fraction = 0.6,
#'     n.neighbors = 20,
#'     reference_pool = "labeled"
#'   )
#'   
#'   # Compare results
#'   table(obj$subclass_nn, obj$subclass_nn_labeled, useNA = "ifany")
#' }
LabelByNearestNeighbors <- function(obj, ident, output_col = NULL, reduction = "pca", 
                                    dims = 1:30, fraction = 0.6, n.neighbors = 20, 
                                    reference_pool = "all") {
  
  library(Seurat)
  library(FNN)
  library(dplyr)
  
  # Validate reference_pool parameter
  if (!reference_pool %in% c("all", "labeled")) {
    stop("reference_pool must be either 'all' or 'labeled'")
  }
  
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
  
  # Determine reference cells based on reference_pool parameter
  if (reference_pool == "all") {
    reference_cells <- c(labeled_cells, unlabeled_cells)
  } else {  # reference_pool == "labeled"
    reference_cells <- labeled_cells
  }
  
  # Find k nearest neighbors for each unlabeled cell
  nn <- get.knnx(
    data = embeddings[reference_cells, , drop = FALSE],      # Reference pool
    query = embeddings[unlabeled_cells, , drop = FALSE],     # Unlabeled cells as query
    k = n.neighbors
  )
  
  # For each unlabeled cell, collect labels from neighbors
  neighbor_idents <- vector("list", length = nrow(nn$nn.index))
  
  for (i in 1:nrow(nn$nn.index)) {
    neighbor_indices <- nn$nn.index[i, ]
    neighbor_cells <- reference_cells[neighbor_indices]
    # Use working copy of labels
    neighbor_idents[[i]] <- labels[match(neighbor_cells, rownames(obj@meta.data))]
  }
  
  # Calculate fraction of neighbors with each label
  all_idents <- unique(labels[match(reference_cells, rownames(obj@meta.data))])
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


# Nearest Neighbor Diagnostics Function
# Call this after each integration round to optimize parameters
diagnose_nn_labeling <- function(obj, ident, reduction = "pca", dims = 1:30, 
                                 reference_pool = "labeled", 
                                 method_col = "method",
                                 spatial_value = "Stereo-seq",
                                 ref_value = "snRNA-seq") {
  
  library(ggplot2)
  library(FNN)
  library(patchwork)
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("NN LABELING DIAGNOSTICS\n"))
  cat(sprintf("========================================\n"))
  cat(sprintf("Ident: %s\n", ident))
  cat(sprintf("Reduction: %s (dims %s)\n", reduction, paste(range(dims), collapse="-")))
  cat(sprintf("Reference pool: %s\n", reference_pool))
  cat(sprintf("========================================\n\n"))
  
  # Get labels
  labels <- obj[[ident]][, 1]
  labels[is.na(labels)] <- "None"
  
  # Identify cell types
  spatial_cells <- colnames(obj)[obj[[method_col]][,1] == spatial_value]
  ref_cells <- colnames(obj)[obj[[method_col]][,1] == ref_value]
  labeled_cells <- rownames(obj@meta.data)[labels != "None"]
  unlabeled_cells <- rownames(obj@meta.data)[labels == "None"]
  
  # Determine if class or subclass based on unique labels
  unique_labels <- unique(labels[labels != "None"])
  is_class_level <- all(unique_labels %in% c("glutamatergic", "gabaergic", "nonneuronal"))
  
  cat(sprintf("Spatial cells: %d\n", length(spatial_cells)))
  cat(sprintf("Reference cells: %d\n", length(ref_cells)))
  cat(sprintf("Labeled cells: %d\n", length(labeled_cells)))
  cat(sprintf("Unlabeled cells: %d (query)\n", length(unlabeled_cells)))
  cat(sprintf("Level: %s\n\n", ifelse(is_class_level, "CLASS", "SUBCLASS")))
  
  # Get embeddings
  coords <- Embeddings(obj, reduction = reduction)
  coords <- coords[, dims, drop = FALSE]
  
  # Determine reference pool
  if (reference_pool == "all") {
    reference_cells <- c(labeled_cells, unlabeled_cells)
  } else {
    reference_cells <- labeled_cells
  }
  
  cat(sprintf("Reference pool size: %d cells\n", length(reference_cells)))
  cat(sprintf("Query pool size: %d cells\n\n", length(unlabeled_cells)))
  
  # Set parameter space based on level and reference_pool
  if (is_class_level) {
    if (reference_pool == "labeled") {
      k_values <- c(25, 50, 100, 200, 500)
      test_fractions <- c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1)
    } else {  # "all"
      k_values <- c(50, 100, 200, 500)
      test_fractions <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
    }
  } else {  # subclass level
    if (reference_pool == "labeled") {
      k_values <- c(10, 25, 50, 100, 200)
      test_fractions <- c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1)
    } else {  # "all"
      k_values <- c(50, 80, 100, 150, 200)
      test_fractions <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
    }
  }
  
  cat(sprintf("Testing K values: %s\n", paste(k_values, collapse=", ")))
  cat(sprintf("Testing fractions: %s\n\n", paste(test_fractions, collapse=", ")))
  
  # Compute nearest neighbors for all K values
  cat("Computing nearest neighbors...\n")
  ref_coords <- coords[reference_cells, , drop = FALSE]
  query_coords <- coords[unlabeled_cells, , drop = FALSE]
  
  nn_results <- list()
  for (k in k_values) {
    nn_results[[as.character(k)]] <- get.knnx(ref_coords, query_coords, k = k)
  }
  
  # === CONSENSUS ANALYSIS ===
  cat("\n--- Consensus Analysis ---\n")
  
  consensus_df <- data.frame()
  
  for (k in k_values) {
    nn_idx <- nn_results[[as.character(k)]]$nn.index
    
    # For each query cell, get labels of K nearest neighbors
    consensus <- sapply(1:nrow(nn_idx), function(i) {
      neighbor_cells <- reference_cells[nn_idx[i, ]]
      neighbor_labels <- labels[match(neighbor_cells, rownames(obj@meta.data))]
      neighbor_labels <- neighbor_labels[neighbor_labels != "None"]  # Exclude "None"
      
      if (length(neighbor_labels) == 0) {
        return(c(NA, 0))
      }
      
      label_table <- table(neighbor_labels)
      max_label <- names(which.max(label_table))
      max_fraction <- max(label_table) / length(neighbor_labels)
      
      c(max_label, max_fraction)
    })
    
    consensus_df <- rbind(consensus_df, data.frame(
      k = k,
      mean_max_fraction = mean(as.numeric(consensus[2,]), na.rm = TRUE),
      median_max_fraction = median(as.numeric(consensus[2,]), na.rm = TRUE),
      pct_above_50 = 100 * sum(as.numeric(consensus[2,]) > 0.5, na.rm = TRUE) / ncol(consensus),
      pct_above_75 = 100 * sum(as.numeric(consensus[2,]) > 0.75, na.rm = TRUE) / ncol(consensus)
    ))
  }
  
  print(consensus_df)
  
  # === DETAILED ANALYSIS ===
  optimal_k <- k_values[ceiling(length(k_values)/2)]  # Middle K value
  cat(sprintf("\n--- Detailed Analysis (K=%d) ---\n", optimal_k))
  
  nn_idx <- nn_results[[as.character(optimal_k)]]$nn.index
  nn_dist <- nn_results[[as.character(optimal_k)]]$nn.dist
  
  query_analysis <- data.frame(
    cell_id = unlabeled_cells,
    dominant_label = NA,
    dominant_fraction = NA,
    mean_nn_dist = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(unlabeled_cells)) {
    neighbor_cells <- reference_cells[nn_idx[i, ]]
    neighbor_labels <- labels[match(neighbor_cells, rownames(obj@meta.data))]
    neighbor_labels <- neighbor_labels[neighbor_labels != "None"]
    
    if (length(neighbor_labels) > 0) {
      label_table <- table(neighbor_labels)
      query_analysis$dominant_label[i] <- names(which.max(label_table))
      query_analysis$dominant_fraction[i] <- max(label_table) / length(neighbor_labels)
    }
    
    query_analysis$mean_nn_dist[i] <- mean(nn_dist[i, ])
  }
  
  cat("\nDominant label distribution:\n")
  print(table(query_analysis$dominant_label, useNA = "ifany"))
  
  cat("\nConsensus strength:\n")
  print(summary(query_analysis$dominant_fraction))
  
  cat("\nNN distance:\n")
  print(summary(query_analysis$mean_nn_dist))
  
  # Check for bimodality in distance distribution
  tryCatch({
    library(diptest)
    dip_test <- dip.test(query_analysis$mean_nn_dist)
    cat(sprintf("\nHartigan's Dip Test for Bimodality: p = %.4f\n", dip_test$p.value))
    if (dip_test$p.value < 0.05) {
      cat("  -> Significant bimodality detected (p < 0.05)\n")
    } else {
      cat("  -> No significant bimodality (p >= 0.05)\n")
    }
  }, error = function(e) {
    cat("\nNote: Install 'diptest' package for bimodality testing\n")
  })
  
  # === VISUALIZATIONS ===
  cat("\nGenerating plots...\n")
  
  # 1. NN Distance Distribution
  p_dist <- ggplot(query_analysis, aes(x = mean_nn_dist)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
    geom_density(aes(y = after_stat(count) * 0.05), color = "red", size = 1) +
    labs(title = sprintf("NN Distance Distribution (K=%d, pool=%s)", optimal_k, reference_pool),
         subtitle = "Distance to K nearest reference cells",
         x = "Mean NN Distance",
         y = "Number of Query Cells") +
    theme_minimal()
  
  print(p_dist)
  
  # 2. Distance vs Consensus scatter
  p_dist_cons <- ggplot(query_analysis, aes(x = mean_nn_dist, y = dominant_fraction)) +
    geom_point(alpha = 0.3, size = 1, color = "steelblue") +
    geom_smooth(method = "loess", color = "red", se = TRUE) +
    labs(title = "Integration Quality vs Consensus",
         subtitle = sprintf("K=%d, pool=%s", optimal_k, reference_pool),
         x = "Mean NN Distance",
         y = "Consensus Strength") +
    theme_minimal()
  
  print(p_dist_cons)
  
  # 3. Distance by Dominant Label
  p_dist_label <- ggplot(query_analysis[!is.na(query_analysis$dominant_label), ], 
                         aes(x = dominant_label, y = mean_nn_dist, fill = dominant_label)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    labs(title = "NN Distance by Dominant Label",
         subtitle = sprintf("K=%d, pool=%s", optimal_k, reference_pool),
         x = "Dominant Label",
         y = "Mean NN Distance") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_dist_label)
  
  # 4. Consensus histogram
  threshold_lines <- test_fractions[c(1, ceiling(length(test_fractions)/2), length(test_fractions))]
  
  p1 <- ggplot(query_analysis, aes(x = dominant_fraction)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
    geom_vline(xintercept = threshold_lines[1], color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = threshold_lines[2], color = "orange", linetype = "dashed", size = 1) +
    geom_vline(xintercept = threshold_lines[3], color = "green", linetype = "dashed", size = 1) +
    labs(title = sprintf("Consensus Strength (K=%d, pool=%s)", optimal_k, reference_pool),
         subtitle = "Fraction of neighbors agreeing on dominant label",
         x = "Dominant Fraction",
         y = "Number of Query Cells") +
    theme_minimal() +
    annotate("text", x = threshold_lines[1], y = Inf, 
             label = sprintf("%.2f", threshold_lines[1]), 
             hjust = -0.1, vjust = 1.5, color = "red", size = 3) +
    annotate("text", x = threshold_lines[2], y = Inf, 
             label = sprintf("%.2f", threshold_lines[2]),
             hjust = -0.1, vjust = 3, color = "orange", size = 3) +
    annotate("text", x = threshold_lines[3], y = Inf, 
             label = sprintf("%.2f", threshold_lines[3]),
             hjust = -0.1, vjust = 4.5, color = "green", size = 3)
  
  print(p1)
  
  # 5. Parameter sweep
  cat("\n--- Parameter Sweep ---\n")
  
  param_results <- data.frame()
  for (k in k_values) {
    nn_idx <- nn_results[[as.character(k)]]$nn.index
    
    for (frac in test_fractions) {
      n_pass <- sum(sapply(1:nrow(nn_idx), function(i) {
        neighbor_cells <- reference_cells[nn_idx[i, ]]
        neighbor_labels <- labels[match(neighbor_cells, rownames(obj@meta.data))]
        neighbor_labels <- neighbor_labels[neighbor_labels != "None"]
        
        if (length(neighbor_labels) == 0) return(FALSE)
        
        label_table <- table(neighbor_labels)
        max(label_table) / length(neighbor_labels) >= frac
      }))
      
      param_results <- rbind(param_results, data.frame(
        k = k,
        fraction = frac,
        n_labeled = n_pass,
        pct_labeled = 100 * n_pass / length(unlabeled_cells)
      ))
    }
  }
  
  print(param_results)
  
  p2 <- ggplot(param_results, aes(x = fraction, y = pct_labeled, color = factor(k))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(title = sprintf("Parameter Optimization (pool=%s)", reference_pool),
         subtitle = "Expected % of query cells labeled",
         x = "Fraction Threshold",
         y = "% Cells Labeled",
         color = "K Neighbors") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal()
  
  print(p2)
  
  # === DISTANCE THRESHOLD ANALYSIS ===
  cat("\n--- Distance-Based Filtering Analysis ---\n")
  
  # Test different distance thresholds
  distance_thresholds <- quantile(query_analysis$mean_nn_dist, 
                                  probs = c(0.25, 0.50, 0.75, 0.90, 0.95), 
                                  na.rm = TRUE)
  
  cat("\nDistance quantiles:\n")
  print(distance_thresholds)
  
  dist_filter_results <- data.frame()
  for (dist_thresh in distance_thresholds) {
    filtered_cells <- query_analysis$mean_nn_dist <= dist_thresh
    n_kept <- sum(filtered_cells, na.rm = TRUE)
    pct_kept <- 100 * n_kept / length(unlabeled_cells)
    mean_consensus_kept <- mean(query_analysis$dominant_fraction[filtered_cells], na.rm = TRUE)
    
    dist_filter_results <- rbind(dist_filter_results, data.frame(
      distance_threshold = dist_thresh,
      n_cells_kept = n_kept,
      pct_cells_kept = pct_kept,
      mean_consensus_kept = mean_consensus_kept
    ))
  }
  
  cat("\nIf filtering by distance threshold:\n")
  print(dist_filter_results)
  
  # === RECOMMENDATIONS ===
  cat("\n========================================\n")
  cat("RECOMMENDATIONS\n")
  cat("========================================\n\n")
  
  # Find parameters for 60-80% labeling
  good_params <- param_results[param_results$pct_labeled >= 60 & 
                                 param_results$pct_labeled <= 80, ]
  
  if (nrow(good_params) > 0) {
    cat("For 60-80% labeling success:\n")
    print(good_params)
    
    # Recommend most conservative (highest fraction)
    best <- good_params[which.max(good_params$fraction), ]
    cat(sprintf("\nRECOMMENDED: fraction=%.2f, n.neighbors=%d (%.1f%% labeled)\n",
                best$fraction, best$k, best$pct_labeled))
    cat(sprintf("            reference_pool='%s'\n\n", reference_pool))
  } else {
    cat("No parameters give 60-80% labeling.\n\n")
    
    best <- param_results[which.max(param_results$pct_labeled), ]
    cat(sprintf("Best available: fraction=%.2f, n.neighbors=%d (%.1f%% labeled)\n",
                best$fraction, best$k, best$pct_labeled))
    cat(sprintf("               reference_pool='%s'\n\n", reference_pool))
  }
  
  # Distance-based recommendation
  median_dist <- median(query_analysis$mean_nn_dist, na.rm = TRUE)
  q75_dist <- quantile(query_analysis$mean_nn_dist, 0.75, na.rm = TRUE)
  
  cat("\nDistance-based filtering options:\n")
  cat(sprintf("  Median distance: %.3f (keeps 50%% of cells)\n", median_dist))
  cat(sprintf("  75th percentile: %.3f (keeps 75%% of cells)\n", q75_dist))
  cat("\nConsider filtering out cells with distance > threshold before labeling\n")
  cat("to focus on well-integrated cells only.\n\n")
  
  # Save to global environment
  result_name <- sprintf("nn_diagnostics_%s_%s", ident, reference_pool)
  assign(result_name, list(
    query_analysis = query_analysis,
    param_results = param_results,
    consensus_df = consensus_df,
    dist_filter_results = dist_filter_results,
    recommended = if(nrow(good_params) > 0) best else param_results[which.max(param_results$pct_labeled), ]
  ), envir = .GlobalEnv)
  
  cat(sprintf("Results saved to: %s\n", result_name))
  cat(sprintf("Access recommended params: %s$recommended\n", result_name))
  cat(sprintf("Access distance analysis: %s$dist_filter_results\n\n", result_name))
  
  invisible(list(
    query_analysis = query_analysis,
    param_results = param_results,
    consensus_df = consensus_df,
    dist_filter_results = dist_filter_results
  ))
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
