#' PlotFeatures
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the visualization family.
#' @param obj (auto) parameter
#' @param features (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family visualization
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
Part of the visualization family.
#' @param obj (auto) parameter
#' @param markers (auto) parameter
#' @param subclass.labels (auto) parameter
#' @param ident.labels (auto) parameter
#' @param savepath (auto) parameter
#' @param ens.id (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family visualization
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
Part of the visualization family.
#' @param obj (auto) parameter
#' @param markers (auto) parameter
#' @param subclass.labels (auto) parameter
#' @param ident.labels (auto) parameter
#' @param savepath (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family visualization
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
Part of the visualization family.
#' @param obj (auto) parameter
#' @param markers (auto) parameter
#' @param polygon.coords (auto) parameter
#' @param savepath (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family visualization
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
