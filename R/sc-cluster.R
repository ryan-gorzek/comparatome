#' ClusterWithSCT
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the clustering family.
#' @param obj (auto) parameter
#' @param resolutions (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family clustering
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
ClusterWithSCT <- function(obj, resolutions) {
  
  obj <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
         RunPCA(npcs = 30, verbose = FALSE) %>%
         FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
         RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
  
  for (res in resolutions) {
    obj <- FindClusters(obj, resolution = res, algorithm = 4, method = "igraph")
  }
  
  return(obj)
  
}


#' PlotClusters
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the clustering family.
#' @param obj (auto) parameter
#' @param group.id (auto) parameter
#' @param save.plots (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family clustering
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotClusters <- function(obj, group.id, save.plots = FALSE) {
  
  if (!missing(group.id)) {
    Idents(obj) <- group.id
  }
  obj$active.ident <- obj@active.ident
  dimplot1 <- DimPlot(obj, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  dimplot2 <- DimPlot(obj, reduction = "umap", group.by = "sample", raster = FALSE, shuffle = TRUE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  cluster.sample <- table(obj$sample, obj$active.ident) %>%
    as.data.frame.matrix() %>%
    rownames_to_column(var = "sample")
  cluster.sample[-1] <- lapply(cluster.sample[-1], function(x) x/sum(x))
  cluster.sample <- cluster.sample %>%
    pivot_longer(
      cols = -c("sample"),
      names_to = "cluster",
      values_to = "count"
    )
  cluster.sample$cluster <- factor(cluster.sample$cluster, levels = unique(cluster.sample$cluster))
  barplot1 <- ggplot(cluster.sample, aes(x=cluster, y=count, fill=sample)) +
    geom_bar(stat="identity") +
    theme_minimal()
  dimplot3 <- DimPlot(obj, reduction = "umap", group.by = "predicted_doublets", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  # Summarize doublets by cluster
  df <- obj[[]]
  df$active.ident <- obj@active.ident
  df_summary <- df %>%
    dplyr::group_by(active.ident, predicted_doublets) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(fraction = count / sum(count))
  # Create the stacked bar plot
  barplot2 <- ggplot(df_summary, aes(x = active.ident, y = fraction, fill = predicted_doublets)) +
    geom_bar(stat = "identity") +
    labs(x = "Clusters", y = "Doublet Fraction", fill = "Value") +
    theme_minimal()
  featplot1 <- FeaturePlot(obj, "nFeature_RNA", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  vlnplot1 <- VlnPlot(obj, "nFeature_RNA")
  featplot2 <- FeaturePlot(obj, "nCount_RNA", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  vlnplot2 <- VlnPlot(obj, "nCount_RNA")
  
  print(dimplot1)
  print(dimplot2)
  print(barplot1)
  print(dimplot3)
  print(barplot2)
  print(featplot1)
  print(vlnplot1)
  print(featplot2)
  print(vlnplot2)
  
  if (save.plots == TRUE){
    ggsave("E:/Opossum_Paper/Figure S1/Opossum_All_Clusters.svg", plot = dimplot1)
    ggsave("E:/Opossum_Paper/Figure S1/Opossum_All_Samples.svg", plot = dimplot2)
    ggsave("E:/Opossum_Paper/Figure S1/Opossum_All_Samples_Bar.svg", plot = barplot1)
    ggsave("E:/Opossum_Paper/Figure S1/Opossum_All_Doublets.svg", plot = dimplot3)
    ggsave("E:/Opossum_Paper/Figure S1/Opossum_All_Doublets_Bar.svg", plot = barplot2)
  }
  
}


#' RemoveCellsByUMAP
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the clustering family.
#' @param seurat_obj (auto) parameter
#' @param umap1_range (auto) parameter
#' @param umap2_range (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family clustering
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
RemoveCellsByUMAP <- function(seurat_obj, umap1_range, umap2_range) {
  # Check if UMAP reduction exists
  if (!"umap" %in% names(seurat_obj@reductions)) {
    stop("UMAP reduction not found in the Seurat object. Please ensure UMAP has been run.")
  }
  
  # Extract UMAP coordinates
  umap_coords <<- Embeddings(seurat_obj, reduction = "umap")
  
  # Find cells within the specified UMAP1 and UMAP2 range
  cells_to_remove <<- which(
    umap_coords[, 1] >= umap1_range[1] & umap_coords[, 1] <= umap1_range[2] &
      umap_coords[, 2] >= umap2_range[1] & umap_coords[, 2] <= umap2_range[2]
  )
  
  # Get cell names to remove
  cells_to_remove <<- rownames(umap_coords)[cells_to_remove]
  
  # Remove cells from Seurat object
  cell_names <- colnames(seurat_obj)
  seurat_obj <- subset(seurat_obj, cells = cell_names[!cell_names %in% cells_to_remove])
  
  return(seurat_obj)
}
