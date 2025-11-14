#' IntegrateObjects
#'
#' Integrate two Seurat objects using SCTransform-based integration workflow.
#' 
#' This function performs cross-species or cross-condition integration of two Seurat objects
#' using Seurat's SCTransform v2 workflow. The integration involves finding shared variation
#' between datasets while preserving biological differences. The function can optionally
#' subsample cells to balance dataset sizes and performs clustering at specified resolutions.
#' 
#' @param seurat_obj1 First Seurat object to integrate.
#' @param seurat_obj2 Second Seurat object to integrate.
#' @param resolutions Numeric vector of clustering resolutions to apply. Default is 1.
#'   Multiple resolutions will create separate clustering columns in metadata.
#' @param nfeatures Integer specifying the number of variable features to use for
#'   integration. Default is 3000. Higher values capture more biological variation but
#'   may include more noise.
#' @param subsample Logical or numeric controlling cell subsampling. If FALSE (default),
#'   no subsampling is performed. If TRUE, the larger dataset is randomly subsampled to
#'   match the smaller dataset size. If numeric, both datasets are subsampled to this
#'   exact number of cells.
#'   
#' @details
#' The integration workflow consists of several steps:
#' 
#' 1. **Feature alignment**: Identifies common genes between the two objects
#'    and subsets both to the shared space.
#'    
#' 2. **Optional subsampling**: If requested, balances cell numbers between datasets
#'    to prevent dominance by larger datasets in the integrated space.
#'    
#' 3. **SCTransform normalization**: Applies SCTransform v2 to each object independently,
#'    which performs variance stabilization and normalization while accounting for
#'    technical variation.
#'    
#' 4. **PCA and UMAP**: Computes initial dimensionality reduction for each object using
#'    30 principal components.
#'    
#' 5. **Integration**: Uses canonical correlation analysis (CCA) to find integration
#'    anchors, then integrates the data while preserving all features. The integration
#'    creates a new "integrated" assay containing batch-corrected expression values.
#'    
#' 6. **Joint analysis**: Performs PCA and UMAP on the integrated data, builds a
#'    shared nearest neighbor graph, and performs Leiden clustering at specified
#'    resolutions.
#' 
#' The resulting integrated object contains:
#' - Original data in the "RNA" assay
#' - SCTransform-normalized data in the "SCT" assay  
#' - Integrated, batch-corrected data in the "integrated" assay
#' - Dimensionality reductions based on the integrated space
#' - Clustering results at specified resolutions
#' 
#' @return Integrated Seurat object containing both datasets in a unified analysis space.
#'   The object includes three assays (RNA, SCT, integrated), PCA and UMAP reductions,
#'   and cluster assignments at each specified resolution.
#' 
#' @export
#' @family integration
#' 
#' @examples
#' \dontrun{
#'  # Basic integration of two species datasets
#'  obj_integrated <- IntegrateObjects(obj_mouse, obj_opossum)
#'  
#'  # Integration with subsampling to balance dataset sizes
#'  obj_integrated <- IntegrateObjects(obj_mouse, obj_opossum, subsample = TRUE)
#'  
#'  # Integration with specific cell count and multiple resolutions
#'  obj_integrated <- IntegrateObjects(
#'    obj_mouse, obj_opossum, 
#'    subsample = 5000,
#'    resolutions = c(0.5, 1, 2),
#'    nfeatures = 5000
#'  )
#'  
#'  # Visualize integration
#'  DimPlot(obj_integrated, group.by = "species")
#'  DimPlot(obj_integrated, group.by = "integrated_snn_res.1")
#' }
IntegrateObjects <- function(seurat_obj1, seurat_obj2, resolutions = 1, nfeatures = 3000, subsample = FALSE) {
  
  # Ensure the two objects have the same set of features (genes)
  common_features <- intersect(rownames(seurat_obj1), rownames(seurat_obj2))
  seurat_obj1 <- subset(seurat_obj1, features = common_features)
  seurat_obj2 <- subset(seurat_obj2, features = common_features)
  
  # If subsample is numeric, subsample both objects to that exact size
  if (!is.logical(subsample)) {
    seurat_obj1 <- seurat_obj1[, sample(colnames(seurat_obj1), subsample)]
    seurat_obj2 <- seurat_obj2[, sample(colnames(seurat_obj2), subsample)]
  } 
  # If subsample is TRUE, balance by subsampling larger to match smaller
  else if (subsample) {
    if (ncol(seurat_obj1) > ncol(seurat_obj2)) {
      seurat_obj1 <- seurat_obj1[, sample(colnames(seurat_obj1), ncol(seurat_obj2))]
    } else if (ncol(seurat_obj2) > ncol(seurat_obj1)) {
      seurat_obj2 <- seurat_obj2[, sample(colnames(seurat_obj2), ncol(seurat_obj1))]
    }
  }
  
  # List of objects to integrate
  seurat_list <- list(seurat_obj1, seurat_obj2)
  
  # Perform SCTransform v2 on each object independently
  seurat_list <- lapply(seurat_list, function(x) {
    x <- SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
      RunPCA(npcs = 30, verbose = FALSE) %>%
      RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
    return(x)
  })
  
  # Select integration features across datasets
  integration_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures)
  
  # Prepare objects for SCT integration
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = integration_features)
  
  # Find integration anchors using CCA
  anchors <- FindIntegrationAnchors(
    object.list = seurat_list, 
    normalization.method = "SCT", 
    anchor.features = integration_features
  )
  
  # Integrate data, preserving all features
  integrated_seurat <- IntegrateData(
    anchorset = anchors, 
    normalization.method = "SCT", 
    features.to.integrate = common_features
  )
  
  # Run PCA and UMAP on the integrated object
  integrated_seurat <- RunPCA(integrated_seurat, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE)
  
  # Perform clustering at specified resolutions
  for (res in resolutions) {
    integrated_seurat <- FindClusters(
      integrated_seurat, 
      resolution = res, 
      algorithm = 4,  # Leiden algorithm
      method = "igraph"
    )
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
