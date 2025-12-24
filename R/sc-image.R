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


#' ImageFeaturePlotPercentile
#'
#' Plot gene expression on spatial coordinates with percentile-based normalization.
#' Expression values are converted to percentile ranks, providing robust normalization
#' that handles outliers and enables cross-species/cross-sample comparisons.
#'
#' @param obj Seurat object with spatial coordinates (FOV image data)
#' @param feature Gene name to plot (must exist in rownames(obj))
#' @param pct.floor Numeric. Percentile floor for color scale (default 0). Values below 
#'   this percentile are clipped to the floor color.
#' @param pct.ceiling Numeric. Percentile ceiling for color scale (default 99). Values 
#'   above this percentile are clipped to the ceiling color. Using 99 instead of 100
#'   prevents outliers from compressing the color scale.
#' @param size Numeric. Point size for plotting (default 2)
#' @param fov Character. FOV name for GetTissueCoordinates (default "COL"). Must match
#'   an image slot in the Seurat object.
#' @param low.color Character. Hex color for low expression (default "#fcf0c2", cream)
#' @param high.color Character. Hex color for high expression (default "#a1020a", dark red)
#' @param show.legend Logical. Whether to display the color legend (default TRUE)
#'
#' @return ggplot object with spatial expression plot
#'
#' @details
#' Percentile normalization provides several advantages over raw expression plotting:
#' \itemize{
#'   \item Robust to outliers: extreme values don't compress the color scale
#'   \item Comparable across datasets: 90th percentile means the same thing regardless
#'     of absolute expression levels
#'   \item Handles zero-inflation: zeros are mapped to the floor percentile, and 
#'     percentile ranks are calculated only among expressing cells. This prevents
#'     the common issue where sparse expression causes all non-zero cells to appear
#'     at high percentiles.
#' }
#'
#' The function:
#' \enumerate{
#'   \item Extracts spatial coordinates using GetTissueCoordinates
#'   \item Fetches expression values using FetchData
#'   \item Separates zero and non-zero expression values
#'   \item Calculates percentile ranks only among expressing cells using ecdf()
#'   \item Maps zeros to floor, applies ceiling to expressing cells
#'   \item Orders points so high-expression cells plot on top
#'   \item Generates plot with coord_fixed(ratio = 1) for square output
#' }
#'
#' Note: Coordinates are plotted with x and y swapped to match typical cortical
#' column orientation (pia at top, white matter at bottom).
#'
#' @export
#' @family spatial-plotting
#'
#' @examples
#' \dontrun{
#'   # Basic usage with defaults
#'   p <- ImageFeaturePlotPercentile(obj.spatial, "Cux2")
#'   
#'   # Adjust ceiling to highlight moderate expression
#'   p <- ImageFeaturePlotPercentile(obj.spatial, "Rorb", pct.ceiling = 95)
#'   
#'   # Custom color scheme (blue-yellow)
#'   p <- ImageFeaturePlotPercentile(
#'     obj.spatial, "Fezf2",
#'     low.color = "#2166AC",
#'     high.color = "#FFEDA0"
#'   )
#'   
#'   # Compare same gene across species with matched normalization
#'   p.mouse <- ImageFeaturePlotPercentile(obj.mouse, "Sox5", pct.ceiling = 99)
#'   p.opossum <- ImageFeaturePlotPercentile(obj.opossum, "Sox5", pct.ceiling = 99)
#' }
ImageFeaturePlotPercentile <- function(obj,
                                       feature,
                                       pct.floor = 0,
                                       pct.ceiling = 99,
                                       size = 2,
                                       fov = "COL",
                                       low.color = "#fcf0c2",
                                       high.color = "#a1020a",
                                       show.legend = TRUE) {
  
  # Validate feature exists
  if (!feature %in% rownames(obj)) {
    stop(paste0("Feature '", feature, "' not found in object"))
  }
  
  # Extract spatial coordinates
  coords <- GetTissueCoordinates(obj, image = fov)
  
  # Extract expression values
  expr <- FetchData(obj, vars = feature)
  
  # Build dataframe
  df <- data.frame(
    x = coords$x,
    y = coords$y,
    expression = expr[[1]]
  )
  
  # Calculate percentiles only among expressing cells to handle zero-inflation
  # Zeros are mapped to floor; non-zeros scaled across remaining range
  expressing <- df$expression > 0
  n_expressing <- sum(expressing)
  
  if (n_expressing > 0) {
    # Calculate percentile ranks only among expressing cells
    expr_nonzero <- df$expression[expressing]
    pct_nonzero <- ecdf(expr_nonzero)(expr_nonzero) * 100
    
    # Initialize: zeros get floor, expressing cells get their percentile
    df$pct <- pct.floor
    df$pct[expressing] <- pct_nonzero
  } else {
    # No expressing cells - all at floor
    df$pct <- pct.floor
  }
  
  # Apply ceiling
  df$pct_scaled <- pmin(df$pct, pct.ceiling)
  
  # Order by expression so high values plot on top
  df <- df[order(df$pct_scaled), ]
  
  # Build plot with swapped coordinates for cortical orientation
  p <- ggplot(df, aes(x = y, y = x, color = pct_scaled)) +
    geom_point(size = size, alpha = 1) +
    scale_color_gradient(
      low = low.color,
      high = high.color,
      limits = c(pct.floor, pct.ceiling),
      name = "Percentile"
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle(feature)
  
  if (!show.legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}


#' ImageFeaturePlotRaw
#'
#' Plot gene expression on spatial coordinates with raw value normalization.
#' Provides direct control over expression cutoffs for visualization.
#'
#' @param obj Seurat object with spatial coordinates (FOV image data)
#' @param feature Gene name to plot (must exist in rownames(obj))
#' @param min.cutoff Numeric. Minimum expression cutoff (default 0). Values below
#'   are clipped to this value.
#' @param max.cutoff Numeric. Maximum expression cutoff (default NULL, uses data max).
#'   Values above are clipped to this value.
#' @param size Numeric. Point size for plotting (default 2)
#' @param fov Character. FOV name for GetTissueCoordinates (default "COL")
#' @param low.color Character. Hex color for low expression (default "#fcf0c2", cream)
#' @param high.color Character. Hex color for high expression (default "#a1020a", dark red)
#' @param show.legend Logical. Whether to display the color legend (default TRUE)
#'
#' @return ggplot object with spatial expression plot
#'
#' @details
#' Unlike ImageFeaturePlotPercentile, this function uses raw expression values
#' for the color scale. This is useful when:
#' \itemize{
#'   \item Comparing absolute expression levels across genes
#'   \item Visualizing data where percentiles would obscure biological signal
#'   \item Creating figures where consistent expression scale is needed
#' }
#'
#' @export
#' @family spatial-plotting
#'
#' @examples
#' \dontrun{
#'   # Plot with automatic max
#'   p <- ImageFeaturePlotRaw(obj.spatial, "Cux2")
#'   
#'   # Clip high expression for better contrast
#'   p <- ImageFeaturePlotRaw(obj.spatial, "Rorb", max.cutoff = 2)
#' }
ImageFeaturePlotRaw <- function(obj,
                                feature,
                                min.cutoff = 0,
                                max.cutoff = NULL,
                                size = 2,
                                fov = "COL",
                                low.color = "#fcf0c2",
                                high.color = "#a1020a",
                                show.legend = TRUE) {
  
  # Validate feature exists
  if (!feature %in% rownames(obj)) {
    stop(paste0("Feature '", feature, "' not found in object"))
  }
  
  # Extract spatial coordinates
  coords <- GetTissueCoordinates(obj, image = fov)
  
  # Extract expression values
  expr <- FetchData(obj, vars = feature)
  
  # Build dataframe
  df <- data.frame(
    x = coords$x,
    y = coords$y,
    expression = expr[[1]]
  )
  
  # Set max cutoff if not provided
  if (is.null(max.cutoff)) {
    max.cutoff <- max(df$expression)
  }
  
  # Apply cutoffs
  df$expr_scaled <- pmax(pmin(df$expression, max.cutoff), min.cutoff)
  
  # Order by expression so high values plot on top
  df <- df[order(df$expr_scaled), ]
  
  # Build plot
  p <- ggplot(df, aes(x = y, y = x, color = expr_scaled)) +
    geom_point(size = size, alpha = 1) +
    scale_color_gradient(
      low = low.color,
      high = high.color,
      limits = c(min.cutoff, max.cutoff),
      name = "Expression"
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle(feature)
  
  if (!show.legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}
