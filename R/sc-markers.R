#' IdentMarkerDict
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the markers family.
#' @param obj (auto) parameter
#' @param subclass.labels (auto) parameter
#' @param ident.labels (auto) parameter
#' @param save.path (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family markers
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the markers family.
#' @param obj (auto) parameter
#' @param subclass.labels (auto) parameter
#' @param save.path (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family markers
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the markers family.
#' @param nested_list (auto) parameter
#' @param subclass (auto) parameter
#' @param clustering_res (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family markers
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the markers family.
#' @param nested_list (auto) parameter
#' @param subclass.col (auto) parameter
#' @param subclass.order (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family markers
#' @examples
#' \dontrun{
#'  # Example usage will be added
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
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the markers family.
#' @param df (auto) parameter
#' @param gene_list (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family markers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
plotGeneFractions <- function(df, gene_list) {
  # Calculate the fraction of genes in each cluster that belong to the gene_list
  fraction_data <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      total_genes = n(),
      matching_genes = sum(gene %in% gene_list),
      fraction = matching_genes / total_genes
    )
  
  # Plot the results
  ggplot(fraction_data, aes(x = cluster, y = fraction)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "",
         x = "Cluster",
         y = "Fraction of Genes") +
    ylim(0, 1) +
    theme_minimal()
}
