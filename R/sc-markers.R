#' IdentMarkerDict
#'
#' Compute differential expression markers for clusters within each subclass across multiple 
#' clustering resolutions. Generates a nested list structure containing marker genes for each 
#' subclass-resolution-cluster combination.
#'
#' @param obj Seurat object with subclass assignments (metadata columns like "subclass.0.2")
#'   and clustering results (columns like "SCT_snn_res.0.2")
#' @param subclass.labels Character vector of subclass names to analyze (e.g., c("IT_A", "Pvalb"))
#' @param ident.labels Character vector of clustering resolution identifiers (e.g., "SCT_snn_res.0.2")
#' @param save.path File path where the marker dictionary RDS will be saved
#'
#' @return Nested list structure: marker.dict[[subclass]][[resolution]] containing 
#'   FindAllMarkers output (data.frame with columns: gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj)
#'
#' @details
#' For each subclass and resolution combination:
#' \enumerate{
#'   \item Subsets object to cells of that subclass
#'   \item Sets identities to the clustering at specified resolution
#'   \item Runs FindAllMarkers with only.pos = TRUE and logfc.threshold = 0.1
#'   \item Stores results in nested list structure
#' }
#'
#' Results are automatically saved to save.path as an RDS file for later retrieval.
#' Skips subclass-resolution pairs with only one cluster.
#'
#' @export
#' @family markers
#'
#' @examples
#' \dontrun{
#'   markers <- IdentMarkerDict(
#'     obj = obj,
#'     subclass.labels = c("IT_A", "IT_B", "L5PT"),
#'     ident.labels = c("SCT_snn_res.0.2", "SCT_snn_res.0.5"),
#'     save.path = "markers/markerdict_clusters.rds"
#'   )
#'   
#'   # Access markers for IT_A at resolution 0.2
#'   it_a_markers <- markers$IT_A$SCT_snn_res.0.2
#' }
IdentMarkerDict <- function(obj, subclass.labels, ident.labels, save.path) {
  
  marker.dict <- list()
  
  for (sbcl in subclass.labels) {
    
    marker.dict[[sbcl]] <- list()
    
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
        if (length(levels(obj.sbcl.id)) > 1) {
          
          marker.dict[[sbcl]][[id]] <- FindAllMarkers(obj.sbcl.id, only.pos = TRUE, logfc.threshold = 0.1)
          
        }
      }
    }
  }
  
  saveRDS(marker.dict, save.path)
  return(marker.dict)
  
}


#' SubclassMarkerDict
#'
#' Compute differential expression markers distinguishing between subclasses.
#' Generates a list containing marker genes that distinguish each subclass from all others.
#'
#' @param obj Seurat object with subclass labels in metadata
#' @param subclass.labels Character, name of metadata column containing subclass assignments
#' @param save.path File path where the marker dictionary RDS will be saved
#'
#' @return List structure: marker.dict[[subclass]] containing FindAllMarkers output 
#'   (data.frame with columns: gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj)
#'   where each "cluster" is actually a subclass
#'
#' @details
#' Sets the specified subclass column as the active identity and runs FindAllMarkers
#' to identify genes enriched in each subclass compared to all others.
#' Uses only.pos = TRUE and logfc.threshold = 0.1.
#'
#' Results saved to save.path as RDS file.
#'
#' @export
#' @family markers
#'
#' @examples
#' \dontrun{
#'   subclass_markers <- SubclassMarkerDict(
#'     obj = obj,
#'     subclass.labels = "subclass",
#'     save.path = "markers/markerdict_subclass.rds"
#'   )
#'   
#'   # View top markers for a subclass
#'   head(subclass_markers$subclass[subclass_markers$subclass$cluster == "IT_A", ])
#' }
SubclassMarkerDict <- function(obj, subclass.labels, save.path) {
  
  marker.dict <- list()
  
  for (sbcl in subclass.labels) {
    
    marker.dict[[sbcl]] <- list()
    
    DefaultAssay(obj) <- "SCT"
    Idents(obj) <- sbcl
      
    marker.dict[[sbcl]] <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.1)

  }
  
  saveRDS(marker.dict, save.path)
  return(marker.dict)
  
}


#' PlotIdentGeneCounts
#'
#' Plot the number of differentially expressed genes per cluster as a function of 
#' log fold-change threshold. Visualizes marker gene specificity and abundance.
#'
#' @param nested_list Marker dictionary from IdentMarkerDict() with structure 
#'   nested_list[[subclass]][[clustering_res]]
#' @param subclass Character, subclass name to plot
#' @param clustering_res Character, clustering resolution to plot (e.g., "SCT_snn_res.0.2")
#'
#' @return ggplot object showing line plots of gene counts vs log2FC threshold for each cluster
#'
#' @details
#' For each cluster within the specified subclass and resolution:
#' \enumerate{
#'   \item Extracts marker genes with their avg_log2FC values
#'   \item Creates a grid of log2FC thresholds from 0.2 to min(max_log2FC, 2)
#'   \item Counts how many genes exceed each threshold
#'   \item Plots as colored lines, one per cluster
#' }
#'
#' Vertical reference lines at 0.2 and 0.5 log2FC highlight common thresholds.
#' Clusters are reverse-sorted for consistent visualization.
#'
#' @export
#' @family markers
#'
#' @examples
#' \dontrun{
#'   markers <- IdentMarkerDict(obj, c("IT_A"), c("SCT_snn_res.0.5"), "markers.rds")
#'   p <- PlotIdentGeneCounts(markers, "IT_A", "SCT_snn_res.0.5")
#'   p + theme(aspect.ratio = 1)
#' }
PlotIdentGeneCounts <- function(nested_list, subclass, clustering_res) {

  # Extract dataframe
  df <- nested_list[[subclass]][[clustering_res]]
  if (!is.null(df)) {
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(df$avg_log2FC, na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 100)
    
    # Calculate the number of genes for each cluster as avg_log2FC varies
    gene_counts <- df %>%
      group_by(cluster) %>%
      do(data.frame(avg_log2FC = log2FC_grid,
                    count = sapply(log2FC_grid, function(x) sum(.$avg_log2FC > x, na.rm = TRUE)))) %>%
      ungroup()
    levels(gene_counts$cluster) <- rev(sort_idents(levels(gene_counts$cluster)))
    # Plot the number of genes
    ggplot(gene_counts, aes(x = avg_log2FC, y = count, color = cluster, group = cluster)) +
      geom_line(size = 1) + # Thicker lines
      labs(title = paste("Number of DE Genes for", subclass, "by", clustering_res),
           x = "avg_log2FC",
           y = "Number of Genes",
           color = "Cluster") +
      theme_minimal() +
      scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
      geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
      annotate("text", x = 0.2, y = max(gene_counts$count) * 0.9, label = "0.2", hjust = -0.2, vjust = -0.5) +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
      annotate("text", x = 0.5, y = max(gene_counts$count) * 0.9, label = "0.5", hjust = -0.2, vjust = -0.5) +
      theme(
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
        axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
      )
  }
}


#' PlotSubclassGeneCounts
#'
#' Plot the number of differentially expressed genes per subclass as a function of 
#' log fold-change threshold. Similar to PlotIdentGeneCounts but for subclass-level markers.
#'
#' @param nested_list Marker dictionary from SubclassMarkerDict() with structure 
#'   nested_list[[subclass.col]]
#' @param subclass.col Character, name of the subclass column used
#' @param subclass.order Character vector specifying display order of subclasses
#'
#' @return ggplot object showing line plots of gene counts vs log2FC threshold for each subclass
#'
#' @details
#' Filters markers to include only those with:
#' \itemize{
#'   \item pct.1 >= 0.2 (expressed in at least 20% of cells in the subclass)
#'   \item p_val_adj < 0.05 (significant after multiple testing correction)
#' }
#'
#' Then generates line plot showing marker counts across log2FC thresholds.
#' Subclasses ordered according to subclass.order parameter.
#'
#' @export
#' @family markers
#'
#' @examples
#' \dontrun{
#'   subclass_markers <- SubclassMarkerDict(obj, "subclass", "markers.rds")
#'   p <- PlotSubclassGeneCounts(
#'     subclass_markers, 
#'     "subclass",
#'     c("IT_A", "IT_B", "L5PT", "L6CT")
#'   )
#'   p + theme(aspect.ratio = 1)
#' }
PlotSubclassGeneCounts <- function(nested_list, subclass.col, subclass.order) {
  
  # Extract dataframe
  df <- nested_list[[subclass.col]] %>%
          filter(pct.1 >= 0.2) %>% 
          filter(p_val_adj < 0.05)
  
  if (!is.null(df)) {
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(df$avg_log2FC, na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 100)
    
    # Calculate the number of genes for each cluster as avg_log2FC varies
    gene_counts <- df %>%
      group_by(cluster) %>%
      do(data.frame(avg_log2FC = log2FC_grid,
                    count = sapply(log2FC_grid, function(x) sum(.$avg_log2FC > x, na.rm = TRUE)))) %>%
      ungroup()
    levels(gene_counts$cluster) <- rev(sort_by_reference(levels(gene_counts$cluster), subclass.order))
    # Plot the number of genes
    ggplot(gene_counts, aes(x = avg_log2FC, y = count, color = cluster, group = cluster)) +
      geom_line(size = 1) + # Thicker lines
      labs(title = paste("Number of DE Genes for Subclasses"),
           x = "avg_log2FC",
           y = "Number of Genes",
           color = "Subclass") +
      theme_minimal() +
      scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
      geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
      annotate("text", x = 0.2, y = max(gene_counts$count) * 0.9, label = "0.2", hjust = -0.2, vjust = -0.5) +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
      annotate("text", x = 0.5, y = max(gene_counts$count) * 0.9, label = "0.5", hjust = -0.2, vjust = -0.5) +
      theme(
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
        axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
      )
  }
}


#' plotGeneFractions
#'
#' Plot the fraction of genes found in a custom gene list among differentially expressed markers.
#'
#' @param df Data.frame of marker genes (from FindAllMarkers)
#' @param gene_list Character vector of genes to search for
#'
#' @return ggplot object showing fraction of markers present in gene_list for each cluster
#'
#' @keywords internal
#' @family markers
#'
#' @examples
#' \dontrun{
#'   canonical_genes <- c("Slc17a6", "Slc17a7", "Sv2b")
#'   p <- plotGeneFractions(markers, canonical_genes)
#' }
plotGeneFractions <- function(df, gene_list) {
  frac_per_cluster <- df %>% 
    group_by(cluster) %>% 
    summarise(fraction = sum(gene %in% gene_list) / n())
  
  levels(frac_per_cluster$cluster) <- rev(sort_idents(levels(frac_per_cluster$cluster)))
  frac_per_cluster$cluster <- factor(frac_per_cluster$cluster, levels = levels(frac_per_cluster$cluster))
  
  ggplot(frac_per_cluster, aes(x = fraction, y = cluster)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Fraction of Genes in Custom List",
         x = "Fraction",
         y = "Cluster")
}
