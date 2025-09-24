#' LabelCells
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the labeling family.
#' @param obj (auto) parameter
#' @param subclass_resolution (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family labeling
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
Part of the labeling family.
#' @param obj (auto) parameter
#' @param ident (auto) parameter
#' @param fraction (auto) parameter
#' @param n.neighbors (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family labeling
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
Part of the labeling family.
#' @param obj (auto) parameter
#' @param y_limits (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family labeling
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
IdentBySample <- function(obj, y_limits = c(0, 0.50)) {
  
  # Assuming your dataframe is named df with columns 'subclass' and 'sample'
  df <- obj[[]]
  df$active.ident <- obj@active.ident
  df$smpl <- df$sample
  specified_order <- levels(obj)

  # Check for NAs and warn the user
  na_samples <- df %>% filter(is.na(active.ident)) %>% count(smpl)
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
    summarize(MedianProportion = median(Proportion))
  
  # Rename the columns for clarity
  colnames(median_proportions) <- c("active.ident", "MedianProportion")
  
  # Convert active.ident to a factor and specify the order
  median_proportions$active.ident <- factor(median_proportions$active.ident, levels = specified_order)
  relative_proportions$active.ident <- factor(relative_proportions$active.ident, levels = specified_order)
  
  # # Set y-axis limits
  # y_limits <- c(0, 0.60) # Modify these values as needed
  
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
#' Auto-generated roxygen skeleton for comparatome.
Part of the labeling family.
#' @param obj (auto) parameter
#' @param subclass.idx (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family labeling
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
