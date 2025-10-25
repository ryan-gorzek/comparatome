#' PlotFeatures
#'
#' Generate UMAP feature plots for a list of genes, printing each plot to the active graphics device.
#' Handles both lists and vectors of gene names, automatically flattening list structures.
#'
#' @param obj Seurat object with UMAP reduction
#' @param features Character vector or named list of genes to plot. If a list, names are ignored 
#'   and all genes are flattened into a single vector.
#'
#' @return NULL (plots printed to graphics device). Silently skips genes not found in the object.
#'
#' @details
#' For each gene in the features:
#' \itemize{
#'   \item Checks if gene exists in rownames(obj)
#'   \item Generates UMAP FeaturePlot with:
#'     - No legend
#'     - Fixed axis limits (-18 to 18 on both axes)
#'     - Equal aspect ratio
#'     - No rasterization
#'   \item Suppresses warnings (typically from missing genes)
#'   \item Prints plot to active device
#' }
#'
#' Useful for quick visualization of marker gene expression patterns.
#'
#' @export
#' @family visualization
#'
#' @examples
#' \dontrun{
#'   # Plot individual genes
#'   PlotFeatures(obj, c("Slc17a6", "Gad1", "Aldh1l1"))
#'   
#'   # Plot from named list (names ignored)
#'   markers <- list(
#'     excitatory = c("Slc17a6", "Slc17a7"),
#'     inhibitory = c("Gad1", "Gad2")
#'   )
#'   PlotFeatures(obj, markers)
#' }
PlotFeatures <- function(obj, features) {
  
  if (class(features) == "list") { features <- unname(unlist(features)) }
  
  for (feat in features) {
    if (feat %in% rownames(obj)) {
      featplot <- FeaturePlot(obj, reduction = "umap", features = feat, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
      featplot <- suppressWarnings(featplot)
      print(featplot)
    }
  }
  
}


#' SaveDotPlots
#'
#' Generate and save dot plots of top marker genes for each subclass and clustering resolution.
#' Creates both fold-change (FC) and percent-difference (PD) ranked dot plots.
#'
#' @param obj Seurat object with subclass assignments and clustering results
#' @param markers Nested marker list from IdentMarkerDict with structure markers[[subclass]][[resolution]]
#' @param subclass.labels Character vector of subclasses to process
#' @param ident.labels Character vector of resolutions to process (e.g., "SCT_snn_res.0.2")
#' @param savepath Base directory for saving plots
#' @param ens.id Ensembl ID prefix for gene name stripping (e.g., "ENSMODG" for opossum)
#'
#' @return NULL (plots saved to disk)
#'
#' @details
#' For each subclass-resolution combination:
#' \enumerate{
#'   \item Subsets object to that subclass
#'   \item For each cluster/type, selects top 20 unique marker genes by:
#'     - avg_log2FC (fold-change) 
#'     - pct.diff (pct.1 - pct.2, percent difference)
#'   \item Generates two dot plots (FC and PD ranked)
#'   \item Strips Ensembl IDs if genes contain the ens.id prefix
#'   \item Saves to: [savepath]/[subclass]/[resolution]/[subclass.col]_DotPlot_FC.png (and _PD.png)
#' }
#'
#' Plots use red-grey color scale without scaling, with vertical gene labels.
#' Only includes markers unique to single clusters (not shared across clusters).
#'
#' @export
#' @family visualization
#'
#' @examples
#' \dontrun{
#'   markers <- IdentMarkerDict(obj, c("IT_A", "Pvalb"), c("SCT_snn_res.0.5"), "markers.rds")
#'   SaveDotPlots(
#'     obj = obj,
#'     markers = markers,
#'     subclass.labels = c("IT_A", "Pvalb"),
#'     ident.labels = c("SCT_snn_res.0.5"),
#'     savepath = "plots/",
#'     ens.id = "ENSMUSG"
#'   )
#' }
SaveDotPlots <- function(obj, markers, subclass.labels, ident.labels, savepath, ens.id) {
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, gsub("/", "", sbcl), "/")
    make_folder(folder_path)
    
    for (id in ident.labels) {
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj) <- sbcl.col
      if (sbcl %in% levels(obj)) {
        
        obj.sbcl.id <- subset(obj, idents = sbcl)
        
        DefaultAssay(obj.sbcl.id) <- "SCT"
        Idents(obj.sbcl.id) <- id
        if (!all(is.na(as.numeric(levels(obj.sbcl.id))))) {
          levels(obj.sbcl.id) <- factor(rev(as.character(sort(as.numeric(levels(obj.sbcl.id))))))
        }
        else { levels(obj.sbcl.id) <- factor(rev(levels(obj.sbcl.id))) }
        if (length(levels(obj.sbcl.id)) > 1) {

          id.path <- paste0(folder_path, gsub("/", "", id), "/")
          make_folder(id.path)
          
          marker_sets <- list(c(1:20))
          for (set in marker_sets) {
            
            all.markers.FC <- list()
            all.markers.PD <- list()
            
            for (type in levels(obj.sbcl.id)) {
              
              all.markers <- markers[[sbcl]][[id]]
              gene.counts <- table(all.markers$gene)
              unique.markers <- all.markers[all.markers$gene %in% names(gene.counts)[gene.counts == 1],]
              type.markers <- unique.markers[unique.markers$cluster == type,]
              type.markers$pct.diff <- type.markers$pct.1 - type.markers$pct.2
              all.markers.FC[[type]] <- top_genes_desc(type.markers, "avg_log2FC", set)
              all.markers.PD[[type]] <- top_genes_desc(type.markers, "pct.diff", set)
              
            }
            
            # make plots
            plot.FC <- DotPlot(obj.sbcl.id, features = rev(all.markers.FC), cols = c("lightgrey", "red"), scale = FALSE) + 
              theme(axis.text.x = element_text(angle = 90)) + 
              theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))
            plot.FC$data$features.plot <- factor(plot.FC$data$features.plot, levels = levels(plot.FC$data$features.plot), 
                                                 labels = unname(strip_if_contains(unlist(rev(all.markers.FC)), ens.id, paste0(ens.id, "000000"))))
            ggsave(paste0(id.path, sbcl.col, "_DotPlot_FC.png"), plot = plot.FC, width = 24, dpi = 300)
            
            plot.PD <- DotPlot(obj.sbcl.id, features = rev(all.markers.PD), cols = c("lightgrey", "red"), scale = FALSE) +
              theme(axis.text.x = element_text(angle = 90)) +
              theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))
            plot.PD$data$features.plot <- factor(plot.PD$data$features.plot, levels = levels(plot.PD$data$features.plot), 
                                                 labels = unname(strip_if_contains(unlist(rev(all.markers.PD)), ens.id, paste0(ens.id, "000000"))))
            ggsave(paste0(id.path, sbcl.col, "_DotPlot_PD.png"), plot = plot.PD, width = 24, dpi = 300)
            
          }
        }
      }
    }
  }
  
}


#' SaveFeaturePlots
#'
#' Generate and save UMAP feature plots for top marker genes in each subclass and resolution.
#' Creates grid layouts of feature plots (5 columns) for both FC and PD ranked genes.
#'
#' @param obj Seurat object with UMAP reduction
#' @param markers Nested marker list from IdentMarkerDict
#' @param subclass.labels Character vector of subclasses to process
#' @param ident.labels Character vector of resolutions to process
#' @param savepath Base directory for saving plots
#'
#' @return NULL (plots saved to disk)
#'
#' @details
#' For each subclass-resolution-type combination:
#' \enumerate{
#'   \item Selects top 20 unique markers by avg_log2FC and pct.diff
#'   \item Creates 5-column grid of FeaturePlots for each gene
#'   \item Applies consistent square axis limits with equal aspect ratio
#'   \item Saves to: [savepath]/[subclass]/[resolution]/[subclass.col]_FeaturePlot_FC_id-[type].png
#' }
#'
#' Handles cases with fewer than 20 genes gracefully.
#' All plots use red-grey color scale and equal-aspect square plotting area.
#'
#' @export
#' @family visualization
#'
#' @examples
#' \dontrun{
#'   markers <- IdentMarkerDict(obj, c("IT_A"), c("SCT_snn_res.0.5"), "markers.rds")
#'   SaveFeaturePlots(
#'     obj = obj,
#'     markers = markers,
#'     subclass.labels = c("IT_A"),
#'     ident.labels = c("SCT_snn_res.0.5"),
#'     savepath = "plots/"
#'   )
#' }
SaveFeaturePlots <- function(obj, markers, subclass.labels, ident.labels, savepath) {
  
  get_plot_limits <- function(plot) {
    gg_build <- ggplot_build(plot)
    xlims <- gg_build$layout$panel_scales_x[[1]]$range$range
    ylims <- gg_build$layout$panel_scales_y[[1]]$range$range
    return(list(xlims = xlims, ylims = ylims))
  }
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, gsub("/", "", sbcl), "/")
    make_folder(folder_path)
    
    for (id in ident.labels) {
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj) <- sbcl.col
      if (sbcl %in% levels(obj)) {
        
        obj.sbcl.id <- subset(obj, idents = sbcl)
        
        DefaultAssay(obj.sbcl.id) <- "SCT"
        Idents(obj.sbcl.id) <- id
        if (!all(is.na(as.numeric(levels(obj.sbcl.id))))) {
          levels(obj.sbcl.id) <- factor(rev(as.character(sort(as.numeric(levels(obj.sbcl.id))))))
        }
        else { levels(obj.sbcl.id) <- factor(rev(levels(obj.sbcl.id))) }
        if (length(levels(obj.sbcl.id)) > 1) {
          
          id.path <- paste0(folder_path, gsub("/", "", id), "/")
          make_folder(id.path)
          
          marker_sets <- list(c(1:20))
          for (set in marker_sets) {
            
            all.markers.FC <- list()
            all.markers.PD <- list()
            
            for (type in levels(obj.sbcl.id)) {
              
              all.markers <- markers[[sbcl]][[id]]
              gene.counts <- table(all.markers$gene)
              unique.markers <- all.markers[all.markers$gene %in% names(gene.counts)[gene.counts == 1],]
              type.markers <- unique.markers[unique.markers$cluster == type,]
              type.markers$pct.diff <- type.markers$pct.1 - type.markers$pct.2
              all.markers.FC[[type]] <- top_genes_desc(type.markers, "avg_log2FC", set)
              all.markers.PD[[type]] <- top_genes_desc(type.markers, "pct.diff", set)
              
              # make feature plots
              features.FC <- unlist(all.markers.FC[[type]])
              features.PD <- unlist(all.markers.PD[[type]])
              
              for (feature_set in list(features.FC, features.PD)) {
                feature_subset <- feature_set[1:min(20, length(feature_set))]
                feature_subset <- feature_subset[!is.na(feature_subset)]
                if (length(feature_subset) > 0) {
                  plots <- FeaturePlot(obj.sbcl.id, features = feature_subset, cols = c("lightgrey", "red"), ncol = 5)
                  
                  # Extract plot limits
                  plot_limits <- get_plot_limits(plots[[1]])
                  xlims <- plot_limits$xlims
                  ylims <- plot_limits$ylims
                  
                  x.range <- diff(xlims)
                  y.range <- diff(ylims)
                  
                  max.range <- max(x.range, y.range)
                  
                  if (x.range < max.range) {
                    xlims <- mean(xlims) + c(-1, 1) * (max.range / 2)
                  }
                  
                  if (y.range < max.range) {
                    ylims <- mean(ylims) + c(-1, 1) * (max.range / 2)
                  }
                  
                  # Adjust the plots with new limits and coord_equal
                  for (i in 1:length(plots$patches$plots)) {
                    plots[[i]] <- plots[[i]] + coord_equal(xlim = xlims, ylim = ylims)
                  }
                  
                  file_prefix <- ifelse(identical(feature_set, features.FC), "FeaturePlot_FC", "FeaturePlot_PD")
                  ggsave(paste0(id.path, sbcl.col, "_", file_prefix, "_id-", gsub("/", "", type), ".png"), plot = plots, width = 24, height = 16, dpi = 300)
                }
              }
            }
          }
        }
      }
    }
  }
}


#' SaveFeaturePlotsByKMeans
#'
#' Generate and save UMAP feature plots for marker genes grouped by k-means clustering,
#' with polygon overlays showing spatial regions of interest.
#'
#' @param obj Seurat object with UMAP reduction
#' @param markers Data.frame containing marker genes with kmeans_cluster column
#' @param polygon.coords Data.frame with columns X..X and Y defining polygon vertices for overlay
#' @param savepath Base path for saving plots (cluster ID will be appended)
#'
#' @return NULL (plots saved to disk)
#'
#' @details
#' Similar to SaveFeaturePlots but groups markers by k-means cluster rather than
#' Seurat clustering. Adds polygon overlays to each feature plot to highlight
#' spatial regions of interest (useful for spatial transcriptomics or anatomical regions).
#'
#' For each k-means cluster:
#' \enumerate{
#'   \item Selects top 20 unique markers by avg_log2FC and pct.diff
#'   \item Creates 5-column grid of FeaturePlots
#'   \item Overlays polygon as dashed black line with points at vertices
#'   \item Saves to: [savepath]_FeaturePlot_FC_KMeans_[cluster].png
#' }
#'
#' @export
#' @family visualization
#'
#' @examples
#' \dontrun{
#'   # Define region of interest
#'   polygon_coords <- data.frame(
#'     X..X = c(100, 200, 200, 100, 100),
#'     Y = c(100, 100, 200, 200, 100)
#'   )
#'   
#'   SaveFeaturePlotsByKMeans(
#'     obj = spatial_obj,
#'     markers = kmeans_markers,
#'     polygon.coords = polygon_coords,
#'     savepath = "plots/spatial_markers"
#'   )
#' }
SaveFeaturePlotsByKMeans <- function(obj, markers, polygon.coords, savepath) {
  
  get_plot_limits <- function(plot) {
    gg_build <- ggplot_build(plot)
    xlims <- gg_build$layout$panel_scales_x[[1]]$range$range
    ylims <- gg_build$layout$panel_scales_y[[1]]$range$range
    return(list(xlims = xlims, ylims = ylims))
  }
  
  clusters <- sort(unique(markers$kmeans_cluster))

    marker_sets <- list(c(1:20))
    for (set in marker_sets) {
        
      all.markers.FC <- list()
      all.markers.PD <- list()
      
      for (cl in clusters) {
        
        all.markers <- markers
        gene.counts <- table(all.markers$gene)
        unique.markers <- all.markers[all.markers$gene %in% names(gene.counts)[gene.counts == 1],]
        type.markers <- unique.markers[unique.markers$kmeans_cluster == cl,]
        type.markers$pct.diff <- type.markers$pct.1 - type.markers$pct.2
        all.markers.FC[[cl]] <- top_genes_desc(type.markers, "avg_log2FC", set)
        all.markers.PD[[cl]] <- top_genes_desc(type.markers, "pct.diff", set)
        
        # make feature plots
        features.FC <- unlist(all.markers.FC[[cl]])
        features.PD <- unlist(all.markers.PD[[cl]])
        
        for (feature_set in list(features.FC, features.PD)) {
          feature_subset <- feature_set[1:min(20, length(feature_set))]
          feature_subset <- feature_subset[!is.na(feature_subset)]
          if (length(feature_subset) > 0) {
            plots <- FeaturePlot(obj, features = feature_subset, cols = c("lightgrey", "red"), ncol = 5)
            
            # Extract plot limits
            plot_limits <- get_plot_limits(plots[[1]])
            xlims <- plot_limits$xlims
            ylims <- plot_limits$ylims
            
            x.range <- diff(xlims)
            y.range <- diff(ylims)
            
            max.range <- max(x.range, y.range)
            
            if (x.range < max.range) {
              xlims <- mean(xlims) + c(-1, 1) * (max.range / 2)
            }
            
            if (y.range < max.range) {
              ylims <- mean(ylims) + c(-1, 1) * (max.range / 2)
            }
            
            # Adjust the plots with new limits and coord_equal
            for (i in 1:length(plots$patches$plots)) {
              plots[[i]] <- plots[[i]] + coord_equal(xlim = xlims, ylim = ylims) + 
                                         geom_point(data = polygon.coords, aes(x = X..X, y = Y), color = "black", size = 2) + 
                                         geom_path(data = polygon.coords, aes(x = X..X, y = Y), color = "black", linetype = "dashed", size = 1)
            }
            file_prefix <- ifelse(identical(feature_set, features.FC), "FeaturePlot_FC", "FeaturePlot_PD")
            ggsave(paste0(savepath, "_", file_prefix, "_KMeans_", cl, ".png"), plot = plots, width = 24, height = 16, dpi = 300)
          }
        }
    }
  }
}
