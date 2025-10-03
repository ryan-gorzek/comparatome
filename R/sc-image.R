#' PlotImageDimGradient
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the spatial-plotting family.
#' @param obj (auto) parameter
#' @param rgb_csv (auto) parameter
#' @param ratio (auto) parameter
#' @param size (auto) parameter
#' @param yticks (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family spatial-plotting
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotImageDimGradient <- function(obj, rgb_csv, ratio = 1, size = 2, yticks = NA) {
  
  # Load RGB colors from CSV
  rgb_data <- read.csv(rgb_csv, row.names = 1)
  
  # Extract spatial coordinates
  coords <- GetTissueCoordinates(obj)
  
  # Create dataframe
  df <- data.frame(x = coords[,1], y = coords[,2], cells = coords$cell)
  
  # Match colors to correct cells
  df <- merge(df, rgb_data, by.x = "cells", by.y = "row.names", all.x = TRUE)
  
  # Ensure proper RGB format
  df$color <- rgb(df$R, df$G, df$B, maxColorValue = 1)
  
  # Reorder points to plot lower intensity ones first
  df <- df[order(df$R + df$G + df$B, decreasing = FALSE), ]
  
  # Generate the plot
  p <- ggplot(df, aes(x = y, y = x, color = color)) + # , color = color
    geom_point(size = size, alpha = 1) +
    guides(color = "none") + 
    scale_color_identity() +  # Directly use RGB colors
    theme_minimal() +
    coord_fixed(ratio = ratio) +
    ylim(yticks[1], yticks[length(yticks)]) +
    scale_y_continuous(breaks = yticks) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  
  return(p)
  
}


#' PlotImageDimClusters
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the spatial-plotting family.
#' @param seurat_obj (auto) parameter
#' @param ident (auto) parameter
#' @param colors_list (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family spatial-plotting
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotImageDimClusters <- function(seurat_obj, ident = "seurat_clusters", colors_list = NA) {
  
  Idents(seurat_obj) <- ident  # Ensure correct identities
  all_clusters <- levels(seurat_obj)
  
  for (cluster in all_clusters) {
    # Assign grey to all clusters except the one being plotted
    cluster_colors <- rep("grey80", length(all_clusters))
    names(cluster_colors) <- all_clusters
    if (class(colors_list) == "list") {
      cluster_colors[cluster] <- colors_list[[cluster]]
    } else {
      cluster_colors[cluster] <- scales::hue_pal()(length(all_clusters))[which(all_clusters == cluster)]
    }
    
    # Generate and print the plot
    plot <- ImageDimPlot(seurat_obj, group.by = ident, 
                         cols = cluster_colors, 
                         size = 3) + ggtitle(paste(cluster))
    print(plot)
  }
  
}


#' MapGenes
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the spatial-plotting family.
#' @param obj (auto) parameter
#' @param mapping_path (auto) parameter
#' @param use_ids (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family spatial-plotting
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
MapGenes <- function(obj, mapping_path, use_ids = FALSE) {
  
  genes.mapping <- read.csv(mapping_path)
  para.idx <- genes.mapping$Gene.stable.ID %in% unique(genes.mapping$Gene.stable.ID[duplicated(genes.mapping$Gene.stable.ID)])
  genes.mapping <- genes.mapping[!para.idx,]
  genes.mapping.other <- as.list(genes.mapping[, 3])
  genes.mapping.self <- as.list(genes.mapping[, 1])
  ids.mapping.self <- as.list(genes.mapping[, 2])
  genes.self <- rownames(obj)
  for (gene in genes.mapping.other) {
    
    idx.other <- which(genes.mapping.other %in% gene)
    
    if (length(idx.other) == 1) {
      
      gene.self <- genes.mapping.self[idx.other]
      id.self <- ids.mapping.self[idx.other]
      
      if ((gene.self == "") | (use_ids == TRUE)) {
        
        idx.self <- which(genes.self %in% id.self)
        genes.self[idx.self] <- gene
        
      } else {
        
        idx.self <- which(genes.self %in% gene.self)
        genes.self[idx.self] <- gene
        
      }
    }
  }
  # Rebuild Seurat object.
  print("Rebuilding object...")
  obj.df <- as.data.frame(as.matrix(obj[["RNA"]]@counts))
  rownames(obj.df) <- genes.self
  obj.temp <- CreateSeuratObject(counts = obj.df, meta.data = obj[[]])
  obj <- obj.temp
  
  return(obj)
  
}
