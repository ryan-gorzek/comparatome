#' PCAProject
#'
#' Project query cells into the PCA space defined by a reference Seurat object.
#' Uses the reference PC loadings to compute PC coordinates for query cells,
#' enabling direct comparison of cells in a shared low-dimensional space without
#' batch correction or integration.
#'
#' @param seurat_query Seurat object containing query cells to project. Must have
#'   scaled data available (run \code{ScaleData} or \code{NormalizeAndPCA} first).
#' @param seurat_reference Seurat object defining the reference PCA space. Must
#'   have PCA computed with \code{RunPCA} or \code{NormalizeAndPCA}.
#' @param reduction Character string specifying the reduction to use from the
#'   reference object (default: "pca"). The corresponding loadings are used for
#'   projection.
#'
#' @return Query Seurat object with the specified reduction replaced by projected
#'   coordinates. The projection uses the reference PC loadings applied to query
#'   scaled expression values.
#'
#' @details
#' **Projection method:**
#' \enumerate{
#'   \item Extracts feature loadings from reference reduction
#'   \item Identifies genes present in both reference loadings and query scaled data
#'   \item Computes projection: \code{query_PCs = t(query_scaled_data) \%*\% reference_loadings}
#'   \item Stores result as a DimReduc object in the query
#' }
#'
#' **Requirements:**
#' \itemize{
#'   \item Both objects must share a substantial set of genes
#'   \item Query object must have scaled data (\code{scale.data} slot populated)
#'   \item Reference object must have the specified reduction with loadings
#'   \item For cross-species work, use shared orthologous genes in both objects
#' }
#'
#' **Use cases:**
#' \itemize{
#'   \item \strong{Cross-species comparison}: Project one species into another's
#'     PC space to visualize evolutionary conservation of transcriptomic structure
#'   \item \strong{Reference mapping}: Project new samples onto a well-characterized
#'     reference atlas
#'   \item \strong{Archetype analysis}: Compare geometric arrangements (e.g.,
#'     tetrahedra from ParTI) across species in a common coordinate system
#'   \item \strong{Batch assessment}: Project held-out samples to evaluate whether
#'     they occupy expected regions of PC space
#' }
#'
#' **Cross-species workflow:**
#' \preformatted{
#'   # 1. Identify shared genes
#'   shared.genes <- intersect(rownames(obj.mouse), rownames(obj.opossum))
#'
#'   # 2. Normalize both with shared genes as variable features
#'   obj.mouse <- NormalizeAndPCA(obj.mouse, features = shared.genes)
#'   obj.opossum <- NormalizeAndPCA(obj.opossum, features = shared.genes)
#'
#'   # 3. Project opossum into mouse PC space
#'   obj.opossum.proj <- PCAProject(obj.opossum, obj.mouse)
#'
#'   # 4. Visualize both species in mouse reference frame
#'   DimPlot(obj.mouse, reduction = "pca", dims = c(1, 2))
#'   DimPlot(obj.opossum.proj, reduction = "pca", dims = c(1, 2))
#' }
#'
#' **Important considerations:**
#' \itemize{
#'   \item Projection is asymmetric: projecting A???B differs from B???A
#'   \item The reference defines the coordinate system; query cells are placed
#'     within it but don't influence the PC directions
#'   \item Large differences in gene expression distributions between query and
#'     reference may lead to projection artifacts
#' }
#'
#' @export
#' @family projection
#'
#' @examples
#' \dontrun{
#'   # Cross-species projection
#'   shared.genes <- intersect(rownames(obj.mouse), rownames(obj.opossum))
#'   obj.mouse <- NormalizeAndPCA(obj.mouse, features = shared.genes)
#'   obj.opossum <- NormalizeAndPCA(obj.opossum, features = shared.genes)
#'
#'   # Project opossum cells into mouse PC space
#'   obj.opossum.in.mouse <- PCAProject(obj.opossum, obj.mouse)
#'
#'   # Compare native vs projected coordinates
#'   DimPlot(obj.opossum, reduction = "pca", group.by = "subclass")
#'   DimPlot(obj.opossum.in.mouse, reduction = "pca", group.by = "subclass")
#'
#'   # Bidirectional projection for symmetric comparison
#'   obj.mouse.in.opossum <- PCAProject(obj.mouse, obj.opossum)
#'
#'   # Export for external analysis (e.g., MATLAB ParTI)
#'   pc_data <- Embeddings(obj.opossum.in.mouse, "pca")[, 1:10]
#'   write.csv(pc_data, "opossum_in_mouse_pcs.csv")
#' }
PCAProject <- function(seurat_query, seurat_reference, reduction = "pca") {
  # Ensure the reference object has the specified reduction
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
  seurat_query[[reduction]] <- CreateDimReducObject(
    embeddings = as.matrix(query_pcs),
    key = paste0(reduction, "_"),
    assay = DefaultAssay(seurat_query)
  )
  
  return(seurat_query)
}


#' MinDistance
#'
#' Calculate the minimum Euclidean distance from each cell to a set of reference
#' points in PC space. This is used to quantify how far cells are from archetype
#' vertices or polygon boundaries in archetype analysis.
#'
#' @param seurat_obj Seurat object with PCA computed.
#' @param reference_df Data frame containing reference point coordinates. Must have
#'   columns corresponding to PC coordinates (see pc1_colname, pc2_colname).
#' @param pc1_colname Character string specifying the column name for PC1 coordinates
#'   in reference_df (default: "PC_1"). Note: In some polygon exports, this may be
#'   named differently (e.g., "X..X" from CSV header cleaning).
#' @param pc2_colname Character string specifying the column name for PC2 coordinates
#'   in reference_df (default: "PC_2"). Similarly may be named "Y" in some exports.
#'
#' @return Seurat object with new metadata column "min_distance_to_reference" containing
#'   the Euclidean distance from each cell to its nearest reference point.
#'
#' @details
#' **Distance calculation:**
#' For each cell, computes the Euclidean distance to all reference points and stores
#' the minimum. This can be used to:
#' \itemize{
#'   \item Identify cells near archetype vertices
#'   \item Assess archetype quality (cells should be near vertices)
#'   \item Correlate distance with mapping quality scores
#' }
#'
#' **Typical workflow with ParTI archetypes:**
#' \preformatted{
#'   # 1. Fit archetype in MATLAB using ParTI
#'   # 2. Export vertex coordinates (arcOrig)
#'   # 3. Import as polygon in R
#'   polygon <- read.csv("archetype_vertices.csv")
#'   
#'   # 4. Calculate distances
#'   obj <- MinDistance(obj, polygon)
#'   
#'   # 5. Visualize: cells near vertices are specialists
#'   FeaturePlot(obj, "min_distance_to_reference", reduction = "pca")
#' }
#'
#' **Note on coordinate systems:**
#' The reference_df coordinates must be in the same coordinate system as the
#' Seurat PCA embeddings. For cross-species projections, ensure consistent
#' transformations (e.g., sign flips) are applied to both cells and reference points.
#'
#' @export
#' @family projection
#'
#' @examples
#' \dontrun{
#'   # Load archetype polygon from MATLAB ParTI output
#'   polygon <- read.csv("triangle_vertices.csv")
#'   
#'   # Close polygon for visualization
#'   polygon <- rbind(polygon, polygon[1, ])
#'   
#'   # Calculate min distance to nearest vertex
#'   obj.CGE <- MinDistance(obj.CGE, polygon, 
#'                          pc1_colname = "X..X", 
#'                          pc2_colname = "Y")
#'   
#'   # Visualize with polygon overlay
#'   DimPlot(obj.CGE, reduction = "pca") +
#'     geom_point(data = polygon, aes(x = X..X, y = Y), size = 3) +
#'     geom_path(data = polygon, aes(x = X..X, y = Y), linetype = "dashed")
#'   
#'   # Correlation: cells near vertices have higher mapping quality?
#'   FeaturePlot(obj.CGE, c("min_distance_to_reference", "predicted.subclass.score"),
#'               reduction = "pca")
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
#' Part of the projection family.
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
