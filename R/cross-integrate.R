#' IntegrateObjects
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the integration family.
#' @param seurat_obj1 (auto) parameter
#' @param seurat_obj2 (auto) parameter
#' @param resolutions (auto) parameter
#' @param nfeatures (auto) parameter
#' @param subsample (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family integration
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
IntegrateObjects <- function(seurat_obj1, seurat_obj2, resolutions = 1, nfeatures = 3000, subsample = FALSE) {
  # Ensure the two objects have the same set of features (genes)
  common_features <- intersect(rownames(seurat_obj1), rownames(seurat_obj2))
  seurat_obj1 <- subset(seurat_obj1, features = common_features)
  seurat_obj2 <- subset(seurat_obj2, features = common_features)
  
  # If subsample is TRUE, randomly subsample the larger object to match the smaller one
  if (!is.logical(subsample)) {
    seurat_obj1 <- seurat_obj1[, sample(colnames(seurat_obj1), subsample)]
    seurat_obj2 <- seurat_obj2[, sample(colnames(seurat_obj2), subsample)]
  } else if (subsample) {
    if (ncol(seurat_obj1) > ncol(seurat_obj2)) {
      seurat_obj1 <- seurat_obj1[, sample(colnames(seurat_obj1), ncol(seurat_obj2))]
    } else if (ncol(seurat_obj2) > ncol(seurat_obj1)) {
      seurat_obj2 <- seurat_obj2[, sample(colnames(seurat_obj2), ncol(seurat_obj1))]
    }
  }
  
  # List of objects to integrate
  seurat_list <- list(seurat_obj1, seurat_obj2)
  
  # Perform SCTransform v2 on each object
    seurat_list <- lapply(seurat_list, function(x) {
      x <- SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
           RunPCA(npcs = 30, verbose = FALSE) %>%
           RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
      return(x)
    })
  # Select integration features
  integration_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures)
  
  # Prepare objects for integration
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = integration_features)
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = integration_features)
  
  # Integrate data
  integrated_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = common_features)
  
  # Run PCA and UMAP on the integrated object
  integrated_seurat <- RunPCA(integrated_seurat, verbose = FALSE) %>%
                       RunUMAP(dims = 1:30, verbose = FALSE) %>%
                       FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE)
  
  for (res in resolutions) {
    integrated_seurat <- FindClusters(integrated_seurat, resolution = res, algorithm = 4, method = "igraph")
  }
  
  return(integrated_seurat)
}


#' PlotIntegration
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the integration family.
#' @param obj (auto) parameter
#' @param integvar (auto) parameter
#' @param integclust (auto) parameter
#' @param subclass.order (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family integration
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotIntegration <- function(obj, integvar, integclust, subclass.order) {
  
  dimplot1 <- DimPlot(obj, reduction = "umap", group.by = integvar, label = FALSE, raster = FALSE, shuffle = TRUE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  print(dimplot1)
  
  for (icl in integclust) {
    
    dimplot1 <- DimPlot(obj, reduction = "umap", split.by = integvar, group.by = icl, label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
    dimplot2 <- DimPlot(obj, reduction = "umap", split.by = integvar, group.by = "subclass", label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
    heatmap1 <- IntegratedClusterOverlapHeatmap(obj, integvar, "subclass", icl, subclass.order)
    dimplot3 <- DimPlot(obj, reduction = "umap", split.by = integvar, group.by = "type", label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
    heatmap2 <- IntegratedClusterOverlapHeatmap(obj, integvar, "type", icl, subclass.order)
    heatmaps3 <- IntegratedClusterMakeupHeatmap(obj, integvar, "type", icl, subclass.order)
    heatmaps4 <- IdentToIntegratedClusterHeatmap(obj, integvar, "type", icl, subclass.order)
    
    print(dimplot1)
    print(dimplot2)
    print(heatmap1)
    print(dimplot3)
    print(heatmap2)
    grid.arrange(grobs = heatmaps3, ncol = 2)
    grid.arrange(grobs = heatmaps4, ncol = 2)
  
  }
}
