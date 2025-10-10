#' ClusterWithSCT
#'
#' Convenience function for SCTransform, RunPCA, FindNeighbors, RunUMAP, and FindClusters at one or more Leiden resolutions.
#' @param obj Seurat object
#' @param resolutions Leiden clustering resolutions, over which the function iterates and stores in metadata as SCT_snn_res.[resolution].
#' @param nPCs Number of principal components (default = 30) to compute (passed to npcs in RunPCA) and use for FindNeighbors and RunUMAP.
#' @return Seurat object after running SCTransform, RunPCA, FindNeighbors, RunUMAP, and FindClusters.
#' @export
#' @family cluster
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
ClusterWithSCT <- function(obj, resolutions, n_PCs = 30) {
  
  # Run Seurat processing pipeline
  obj <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
         RunPCA(npcs = n_PCs, verbose = FALSE) %>%
         FindNeighbors(reduction = "pca", dims = 1:n_PCs, verbose = FALSE) %>%
         RunUMAP(reduction = "pca", dims = 1:n_PCs, verbose = FALSE)
  # CLuster at each specified resolution, automatically stored as SCT_snn_res.[resolution]
  for (res in resolutions) {
    obj <- FindClusters(obj, resolution = res, algorithm = 4, method = "igraph")
  }
  return(obj)
}


#' PlotClusters
#'
#' Generate useful post-clustering plots for preprocessing.
#' Data are grouped by obj@active.ident, unless group.id (corresponding to a metadata column) is specified.
#' Plots include (default grouping is active.ident or group.id, unless listed):
#' 1. DimPlot (name in output list: dimplot_cluster)
#' 2. DimPlot grouped by 'sample' (dimplot_sample)
#' 3. Stacked bar plot of relative 'sample' frequency for each grouping value (barplot_sample)
#' 4. DimPlot grouped by 'predicted_doublets' from ScrubletR run in PreprocessData (dimplot_doublet)
#' 5. Stacked bar plot of relative 'predicted_doublets' frequency for each grouping value (barplot_doublet)
#' 6. FeaturePlot scaled by 'nFeature_RNA' (featplot_nfeature)
#' 7. VlnPlot of 'nFeature_RNA' for each grouping value (vlnplot_nfeature)
#' 8. FeaturePlot scaled by 'nCount_RNA' (featplot_ncount)
#' 9. VlnPlot of 'nCount_RNA' for each grouping value (vlnplot_ncount)
#' @param obj Seurat object, after running ClusterWithSCT.
#' @param group.id Optional Metadata column set as active ident
#' @param x_lim X-axis limits for DimPlots, default c(-20, 20)
#' @param y_lim Y-axis limits for DimPlots, default c(-20, 20)
#' @return List of plots, named as in description of indivdual plots.
#' @export
#' @family cluster
#' @examples
#' \dontrun{
#'  obj.opossum <- ClusterWithSCT(obj.opossum, 1)
#'  clust.plots <- PlotClusters(obj.opossum)
#'  save.plots <- list("dimplot_cluster" = "S1E_All-Cluster-DimPlot",
#'                     "dimplot_sample" = "S1F_All-Sample-DimPlot",
#'                     "barplot_sample" = "S1G_All-Sample-BarPlot",
#'                     "barplot_doublet" = "S1H_All-Doublet-BarPlot")
#'  for (sp in names(save.plots)) {
#'    SavePNGandSVG(clust.plots[[sp]], dir.list$figS1$plots, save.plots[[sp]])
#'  }
#' }
PlotClusters <- function(obj, group.id, x_lim = c(-20, 20), y_lim = c(-20, 20)) {
  
  # Set Idents with group.id, if specified
  if (!missing(group.id)) {
    Idents(obj) <- group.id
  }
  # Store active ident for grouping
  active.ident <- obj@active.ident
  # Generate plots:
  # dimplot_cluster
  dimplot_cluster <- DimPlot(obj, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend() + xlim(x_lim[1], x_lim[2]) + ylim(y_lim[1], y_lim[2]) + coord_equal()
  # dimplot_sample
  dimplot_sample <- DimPlot(obj, reduction = "umap", group.by = "sample", raster = FALSE, shuffle = TRUE) + xlim(x_lim[1], x_lim[2]) + ylim(y_lim[1], y_lim[2]) + coord_equal()
  # barplot_sample
  cluster.sample <- table(obj$sample, active.ident) %>%
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
  barplot_sample <- ggplot(cluster.sample, aes(x=cluster, y=count, fill=sample)) +
    geom_bar(stat="identity") +
    theme_minimal()
  # dimplot_doublet
  dimplot_doublet <- DimPlot(obj, reduction = "umap", group.by = "predicted_doublets", raster = FALSE) + xlim(x_lim[1], x_lim[2]) + ylim(y_lim[1], y_lim[2]) + coord_equal()
  # barplot_doublet
  df <- obj[[]]
  df$active.ident <- obj@active.ident
  df_summary <- df %>%
    dplyr::group_by(active.ident, predicted_doublets) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(fraction = count / sum(count))
  barplot_doublet <- ggplot(df_summary, aes(x = active.ident, y = fraction, fill = predicted_doublets)) +
    geom_bar(stat = "identity") +
    labs(x = "Clusters", y = "Doublet Fraction", fill = "Value") +
    theme_minimal()
  # nFeature and nCount plots
  featplot_nfeature <- FeaturePlot(obj, "nFeature_RNA", raster = FALSE) + xlim(x_lim[1], x_lim[2]) + ylim(y_lim[1], y_lim[2]) + coord_equal()
  vlnplot_nfeature <- VlnPlot(obj, "nFeature_RNA")
  featplot_ncount <- FeaturePlot(obj, "nCount_RNA", raster = FALSE) + xlim(x_lim[1], x_lim[2]) + ylim(y_lim[1], y_lim[2]) + coord_equal()
  vlnplot_ncount <- VlnPlot(obj, "nCount_RNA")
  # Print plots for notebook visualization
  print(dimplot_cluster)
  print(dimplot_sample)
  print(barplot_sample)
  print(dimplot_doublet)
  print(barplot_doublet)
  print(featplot_nfeature)
  print(vlnplot_nfeature)
  print(featplot_ncount)
  print(vlnplot_ncount)
  # Return plots in list for creating figures
  return(list(dimplot_cluster = dimplot_cluster,
              dimplot_sample = dimplot_sample,
              barplot_sample = barplot_sample,
              dimplot_doublet = dimplot_doublet,
              barplot_doublet = barplot_doublet,
              featplot_nfeature = featplot_nfeature,
              vlnplot_nfeature = vlnplot_nfeature,
              featplot_ncount = featplot_ncount,
              vlnplot_ncount = vlnplot_ncount))
}


#' RemoveCellsByUMAP
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the clustering family.
#' @param seurat_obj (auto) parameter
#' @param umap1_range (auto) parameter
#' @param umap2_range (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family cluster
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
