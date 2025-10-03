#' PlotIdentDEIntersection
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the de-overlap family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param all_genes1 (auto) parameter
#' @param all_genes2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass1 (auto) parameter
#' @param subclass2 (auto) parameter
#' @param ident1 (auto) parameter
#' @param ident2 (auto) parameter
#' @param log2FC_threshold (auto) parameter
#' @param percentage (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-overlap
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotIdentDEIntersection <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass1, subclass2, ident1, ident2, log2FC_threshold, percentage = FALSE) {
  # Extract dataframes
  df1 <- list1[[subclass1]][[ident1]]
  df2 <- list2[[subclass2]][[ident2]]
  
  # Filter by log2FC_threshold
  df1 <- df1 %>% filter(avg_log2FC > log2FC_threshold)
  df2 <- df2 %>% filter(avg_log2FC > log2FC_threshold)
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Create grid
  cluster_combinations <- expand.grid(
    Cluster1 = unique(df1$cluster), 
    Cluster2 = unique(df2$cluster)
  )
  
  # Initialize an empty data frame to store intersection counts or percentages
  intersection_counts <- data.frame()
  
  # Calculate intersections for each cluster pair
  for (i in 1:nrow(cluster_combinations)) {
    cluster1 <- cluster_combinations$Cluster1[i]
    cluster2 <- cluster_combinations$Cluster2[i]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Find intersecting DE genes for the specific cluster pair
    intersecting_de_genes <- intersect(df1_filtered$gene, df2_filtered$gene)
    
    # Count the number of intersecting DE genes
    count <- length(intersecting_de_genes)
    
    if (percentage) {
      # Calculate total intersecting DE genes for percentage calculation
      # total_intersecting_genes <<- length(intersect(df1_filtered$gene, intersecting_genes)) + length(intersect(df2_filtered$gene, intersecting_genes))
      total_intersecting_genes <- length(union(intersect(df1_filtered$gene, intersecting_genes), intersect(df2_filtered$gene, intersecting_genes)))
      # Calculate percentage
      percent <- (count / total_intersecting_genes) * 100
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = percent))
    } else {
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = count))
    }
  }
  
  fill_label <- if (percentage) "Percentage of DE Genes" else "# DE Genes"
  
  intersection_counts$Cluster1 <- factor(intersection_counts$Cluster1, levels = rev(sort_idents(levels(intersection_counts$Cluster1))))
  intersection_counts$Cluster2 <- factor(intersection_counts$Cluster2, levels = rev(sort_idents(levels(intersection_counts$Cluster2))))
  
  # Plot heatmap
  ggplot(intersection_counts, aes(x = Cluster1, y = Cluster2, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Value, 2)), color = "black") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, max(intersection_counts$Value))) +
    scale_y_discrete(limits = rev(levels(factor(intersection_counts$Cluster2)))) + # Reverse y-axis order
    labs(
      title = paste0("Shared DE Genes (log2FC > ", log2FC_threshold, ")"), 
      x = sample.name1,
      y = sample.name2,
      fill = fill_label
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    ) # Make cells square
}


#' SaveSubclassDEIntersectionGenes
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the de-overlap family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param all_genes1 (auto) parameter
#' @param all_genes2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param log2FC_thresholds (auto) parameter
#' @param output_path (auto) parameter
#' @param percentage (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-overlap
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
SaveSubclassDEIntersectionGenes <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, log2FC_thresholds, output_path, percentage = FALSE) {
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
   # Filter dataframes based on thresholds and subclass order
  df1 <- list1$subclass %>%
    filter(pct.1 >= 0.2, p_val_adj < 0.05, gene %in% intersecting_genes) %>%
    filter(cluster %in% subclass.order)
  df1$cluster <- factor(df1$cluster, levels = subclass.order)
  
  df2 <- list2$subclass %>%
    filter(pct.1 >= 0.2, p_val_adj < 0.05, gene %in% intersecting_genes) %>%
    filter(cluster %in% subclass.order)
  df2$cluster <- factor(df2$cluster, levels = subclass.order)
  
  # Create grid of subclass combinations
  cluster_combinations <- expand.grid(
    Cluster1 = unique(df1$cluster), 
    Cluster2 = unique(df2$cluster)
  )
  
  # Iterate over log2FC thresholds
  for (log2FC_threshold in log2FC_thresholds) {
    # Filter by log2FC threshold
    df1_filtered <- df1 %>% filter(avg_log2FC > log2FC_threshold)
    df2_filtered <- df2 %>% filter(avg_log2FC > log2FC_threshold)
    
    # Initialize an empty dataframe for intersection counts or percentages
    intersection_counts <- data.frame()
    
    # Process each cluster pair
    for (i in 1:nrow(cluster_combinations)) {
      cluster1 <- cluster_combinations$Cluster1[i]
      cluster2 <- cluster_combinations$Cluster2[i]
      
      # Filter for specific cluster pair
      df1_cluster <- df1_filtered %>% filter(cluster == cluster1)
      df2_cluster <- df2_filtered %>% filter(cluster == cluster2)
      
      # Find intersecting and unique DE genes
      intersecting_de_genes <- intersect(df1_cluster$gene, df2_cluster$gene)
      de_genes_1 <- setdiff(df1_cluster$gene, intersecting_de_genes)
      de_genes_2 <- setdiff(df2_cluster$gene, intersecting_de_genes)
      
      # Count intersections or calculate percentages
      count <- length(intersecting_de_genes)
      if (percentage) {
        total_genes <- length(union(df1_cluster$gene, df2_cluster$gene))
        percent <- (count / total_genes) * 100
        intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = percent))
      } else {
        intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = count))
      }
      
      # Format log2FC threshold for filename
      threshold_label <- formatC(log2FC_threshold, format = "e", digits = 1)
      threshold_label <- gsub("\\.", "e", threshold_label)
      
      # Write output files
      intersect_file <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_intersecting_genes_", threshold_label, ".txt")
      write.table(intersecting_de_genes, file = intersect_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      de_file_1 <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_", tolower(sample.name1),"_genes_", threshold_label, ".txt")
      write.table(de_genes_1, file = de_file_1, quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      de_file_2 <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_", tolower(sample.name2),"_genes_", threshold_label, ".txt")
      write.table(de_genes_2, file = de_file_2, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
}


#' PlotSubclassDEIntersectionHeatmap
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the de-overlap family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param all_genes1 (auto) parameter
#' @param all_genes2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param log2FC_threshold (auto) parameter
#' @param output_path (auto) parameter
#' @param percentage (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-overlap
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassDEIntersectionHeatmap <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, log2FC_threshold, output_path, percentage = FALSE) {
  # Extract dataframes
  df1 <- list1$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  df2 <- list2$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)

  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Filter by log2FC_threshold
  df1 <- df1 %>%
           filter(gene %in% intersecting_genes) %>%
           filter(avg_log2FC > log2FC_threshold & cluster %in% subclass.order)
  subclasses <- subclass.order[subclass.order %in% df1$cluster]
  df1$cluster <- factor(df1$cluster, levels = subclasses)
  df2 <- df2 %>%
           filter(gene %in% intersecting_genes) %>%
           filter(avg_log2FC > log2FC_threshold & cluster %in% subclass.order)
  subclasses <- subclass.order[subclass.order %in% df2$cluster]
  df2$cluster <- factor(df2$cluster, levels = subclasses)

  # Create grid
  cluster_combinations <<- expand.grid(
    Cluster1 = unique(df1$cluster), 
    Cluster2 = unique(df2$cluster)
  )
  
  # Initialize an empty data frame to store intersection counts or percentages
  intersection_counts <- data.frame()
  
  # Calculate intersections for each cluster pair
  for (i in 1:nrow(cluster_combinations)) {
    cluster1 <- cluster_combinations$Cluster1[i]
    cluster2 <- cluster_combinations$Cluster2[i]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Find intersecting DE genes for the specific cluster pair
    intersecting_de_genes <- intersect(df1_filtered$gene, df2_filtered$gene)
    de_genes_1 <- as.character(df1_filtered$gene[(df1_filtered$gene %in% intersecting_de_genes) == FALSE])
    de_genes_2 <- as.character(df2_filtered$gene[(df2_filtered$gene %in% intersecting_de_genes) == FALSE])
    
    # Count the number of intersecting DE genes
    count <- length(intersecting_de_genes)
    
    if (percentage) {
      # Calculate total intersecting DE genes for percentage calculation
      total_intersecting_genes <- length(union(df1_filtered$gene, df2_filtered$gene))
      # Calculate percentage
      percent <- (count / total_intersecting_genes) * 100
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = percent))
    } else {
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = count))
    }
    
    # Save intersecting genes to a text file
    filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_intersecting_genes.txt")
    write.table(intersecting_de_genes, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
    filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_", tolower(sample.name1), "_genes.txt")
    write.table(de_genes_1, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
    filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_", tolower(sample.name2), "_genes.txt")
    write.table(de_genes_2, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  fill_label <- if (percentage) "Percentage of DE Genes" else "# DE Genes"
  
  intersection_counts$Cluster1 <- factor(intersection_counts$Cluster1, levels = rev(sort_by_reference(levels(intersection_counts$Cluster1), subclass.order)))
  intersection_counts$Cluster2 <- factor(intersection_counts$Cluster2, levels = rev(sort_by_reference(levels(intersection_counts$Cluster2), subclass.order)))
  
  # Plot heatmap
  ggplot(intersection_counts, aes(x = Cluster1, y = Cluster2, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Value, 1)), color = "black") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 20), oob = scales::squish) +
    scale_y_discrete(limits = rev(levels(factor(intersection_counts$Cluster2)))) + # Reverse y-axis order
    labs(
      title = paste0("Shared DE Genes (log2FC > ", log2FC_threshold, ")"), 
      x = sample.name1,
      y = sample.name2,
      fill = fill_label
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = length(levels(intersection_counts$Cluster2)) / length(levels(intersection_counts$Cluster1)),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis tick labels
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    ) # Make cells square
}


#' PlotSubclassDEIntersectionHeatmapTopX
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the de-overlap family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param all_genes1 (auto) parameter
#' @param all_genes2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param top_genes_threshold (auto) parameter
#' @param output_path (auto) parameter
#' @param percentage (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-overlap
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassDEIntersectionHeatmapTopX <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, top_genes_threshold, output_path, percentage = FALSE) {
  # Extract dataframes
  df1 <- list1$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  df2 <- list2$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Filter by top X genes
  df1 <- df1 %>% 
    filter(gene %in% intersecting_genes) %>%
    group_by(cluster) %>%
    top_n(n = top_genes_threshold, wt = avg_log2FC) %>% 
    filter(cluster %in% subclass.order)
  subclasses <- subclass.order[subclass.order %in% df1$cluster]
  df1$cluster <- factor(df1$cluster, levels = subclasses)
  
  df2 <- df2 %>% 
    filter(gene %in% intersecting_genes) %>%
    group_by(cluster) %>%
    top_n(n = top_genes_threshold, wt = avg_log2FC) %>% 
    filter(cluster %in% subclass.order)
  subclasses <- subclass.order[subclass.order %in% df2$cluster]
  df2$cluster <- factor(df2$cluster, levels = subclasses)

  # Create grid
  cluster_combinations <<- expand.grid(
    Cluster1 = unique(df1$cluster), 
    Cluster2 = unique(df2$cluster)
  )
  
  # Initialize an empty data frame to store intersection counts or percentages
  intersection_counts <- data.frame()
  
  # Calculate intersections for each cluster pair
  for (i in 1:nrow(cluster_combinations)) {
    cluster1 <- cluster_combinations$Cluster1[i]
    cluster2 <- cluster_combinations$Cluster2[i]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Find intersecting DE genes for the specific cluster pair
    intersecting_de_genes <- intersect(df1_filtered$gene, df2_filtered$gene)
    
    # Count the number of intersecting DE genes
    count <- length(intersecting_de_genes)
    
    if (percentage) {
      # Calculate total intersecting DE genes for percentage calculation
      total_intersecting_genes <- length(union(df1_filtered$gene, df2_filtered$gene))
      # Calculate percentage
      percent <- (count / total_intersecting_genes) * 100
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = percent))
    } else {
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = count))
    }
    
    # Save intersecting genes to a text file
    filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_intersecting_genes.txt")
    write.table(intersecting_de_genes, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
    # filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_all_genes.txt")
    # write.table(all_de_genes, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  fill_label <- if (percentage) "Percentage of DE Genes" else "# DE Genes"
  
  intersection_counts$Cluster1 <- factor(intersection_counts$Cluster1, levels = rev(sort_by_reference(levels(intersection_counts$Cluster1), subclass.order)))
  intersection_counts$Cluster2 <- factor(intersection_counts$Cluster2, levels = rev(sort_by_reference(levels(intersection_counts$Cluster2), subclass.order)))
  
  # Plot heatmap
  ggplot(intersection_counts, aes(x = Cluster1, y = Cluster2, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Value, 1)), color = "black") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 20), oob = scales::squish) + # max(intersection_counts$Value)
    scale_y_discrete(limits = rev(levels(factor(intersection_counts$Cluster2)))) + # Reverse y-axis order
    labs(
      title = paste0("Shared DE Genes (Top ", top_genes_threshold, " Genes)"), 
      x = sample.name1,
      y = sample.name2,
      fill = fill_label
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = length(levels(intersection_counts$Cluster2)) / length(levels(intersection_counts$Cluster1)),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis tick labels
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    ) # Make cells square
}
