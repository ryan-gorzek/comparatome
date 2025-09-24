#' PCAProject
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the projection family.
#' @param seurat_query (auto) parameter
#' @param seurat_reference (auto) parameter
#' @param reduction (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family projection
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PCAProject <- function(seurat_query, seurat_reference, reduction = "pca") {
  # Ensure the reference object has the specified reduction (e.g., PCA)
  if (!reduction %in% names(seurat_reference@reductions)) {
    stop(paste0("Reduction '", reduction, "' not found in the reference object."))
  }
  
  # Extract the PC loadings from the reference object
  reference_loadings <- seurat_reference[[reduction]]@feature.loadings
  
  # Extract the scaled data from the query object
  query_scaled_data <- GetAssayData(seurat_query, slot = "scale.data")
  
  # Ensure the features (genes) overlap between the query and reference
  common_features <- intersect(rownames(reference_loadings), rownames(query_scaled_data))
  if (length(common_features) < 1) {
    stop("No overlapping features found between the reference and query objects.")
  }
  
  # Subset the loadings and scaled data to the common features
  reference_loadings <- reference_loadings[common_features, , drop = FALSE]
  query_scaled_data <- query_scaled_data[common_features, , drop = FALSE]
  
  # Project the query cells onto the reference PCs
  query_pcs <- t(query_scaled_data) %*% reference_loadings
  
  # Add the projected PCs back to the query object
  seurat_query[[reduction]] <- CreateDimReducObject(embeddings = as.matrix(query_pcs), key = paste0(reduction, "_"), assay = DefaultAssay(seurat_query))
  
  return(seurat_query)
}


#' MinDistance
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the projection family.
#' @param seurat_obj (auto) parameter
#' @param reference_df (auto) parameter
#' @param pc1_colname (auto) parameter
#' @param pc2_colname (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family projection
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
MinDistance <- function(seurat_obj, reference_df, pc1_colname = "PC_1", pc2_colname = "PC_2") {
  # Extract PC1 and PC2 coordinates from the Seurat object
  cell_coords <- as.data.frame(Embeddings(seurat_obj, "pca")[, c(pc1_colname, pc2_colname)])
  
  # Function to calculate the Euclidean distance between two points
  euclidean_dist <- function(x1, y1, x2, y2) {
    sqrt((x1 - x2)^2 + (y1 - y2)^2)
  }
  
  # Initialize a vector to store the minimum distances
  min_distances <- numeric(nrow(cell_coords))
  
  # Calculate the minimum distance for each cell
  for (i in 1:nrow(cell_coords)) {
    distances <- apply(reference_df, 1, function(ref_point) {
      euclidean_dist(cell_coords[i, "PC_1"], cell_coords[i, "PC_2"], ref_point["X..X"], ref_point["Y"])
    })
    min_distances[i] <- min(distances)
  }
  
  # Add the distances to the Seurat object's metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = min_distances, col.name = "min_distance_to_reference")
  
  return(seurat_obj)
}


#' PlotVarianceExplainedByPCs
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the projection family.
#' @param seurat_objects (auto) parameter
#' @param object_names (auto) parameter
#' @param num_pcs (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family projection
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotVarianceExplainedByPCs <- function(seurat_objects, object_names, num_pcs = 10) {
  library(ggplot2)
  library(matrixStats)
  
  # Check if inputs are vectors
  if (!is.list(seurat_objects)) {
    seurat_objects <- list(seurat_objects)
  }
  if (!is.character(object_names)) {
    object_names <- as.character(object_names)
  }
  
  # Initialize a dataframe to hold the variance explained data
  variance_data <- data.frame(PC = integer(), Variance_Explained = numeric(), Species = character())
  
  # Loop through each Seurat object to calculate the variance explained
  for (i in seq_along(seurat_objects)) {
    obj <- seurat_objects[[i]]
    obj_name <- object_names[i]
    
    # Calculate the total variance using rowVars on the scale.data
    total_variance <- sum(rowVars(obj@assays$SCT@scale.data))
    
    # Calculate the variance explained by each PC
    variance_explained <- obj@reductions$pca@stdev^2 / total_variance * 100
    
    # Select the top PCs
    top_variance <- variance_explained[1:num_pcs]
    
    # Create a temporary dataframe for this object
    temp_df <- data.frame(
      PC = 1:num_pcs,
      Variance_Explained = top_variance,
      Species = obj_name
    )
    
    # Add this data to the main dataframe
    variance_data <- rbind(variance_data, temp_df)
  }
  
  # Plot the scree plot
  ggplot(variance_data, aes(x = PC, y = Variance_Explained, color = Species, group = Species)) +
    geom_line() +
    geom_point() +
    labs(title = paste("Scree Plot of Top", num_pcs, "PCs"),
         x = "Principal Component",
         y = "% Variance Explained") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "red")) +
    scale_x_continuous(breaks = 1:num_pcs, labels = 1:num_pcs)
}
