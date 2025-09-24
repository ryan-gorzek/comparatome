#' PlotSubclassDEIntersectionTotalCDF
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the de-curves family.
#' @param list (auto) parameter
#' @param ortho_genes (auto) parameter
#' @param subclass.order (auto) parameter
#' @param subclass_colors (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-curves
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassDEIntersectionTotalCDF <- function(list, ortho_genes, subclass.order, subclass_colors) {
  # Extract dataframes
  df <- list$subclass %>%
          filter(pct.1 >= 0.2) %>% 
          filter(p_val_adj < 0.05)
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  # Calculate cumulative distribution for each specified cluster pair
  for (i in seq_along(subclass.order)) {
    sbcl <- subclass.order[[i]]
    
    # Filter dataframes for the specific cluster pair
    df_filtered <- df %>% filter(cluster == sbcl)
    
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(df_filtered$avg_log2FC, na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 50)
    
    count <- c()
    for (l in log2FC_grid) {
      de_genes <- df_filtered$gene[df_filtered$avg_log2FC > l]
      de_genes_ortho <- de_genes[de_genes %in% ortho_genes]
      count <- c(count, length(de_genes_ortho) / length(de_genes))
    }
    
    # Calculate the cumulative number of shared genes as avg_log2FC varies
    shared_gene_counts <- data.frame(
      avg_log2FC = log2FC_grid,
      count = count,
      Cluster = sbcl,
      Color = subclass_colors[i]
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
  }
  
  cdf_data$Cluster <- factor(cdf_data$Cluster, levels = subclass.order)
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count, color = Cluster, group = Cluster)) +
    geom_line(size = 1) + # Use line size
    scale_color_manual(values = setNames(subclass_colors, subclass.order), guide = guide_legend(reverse = TRUE)) + # Set the colors using subclass_colors
    labs(title = paste0(""),
         x = "avg_log2FC",
         y = "Fraction of DE Genes with a 1:1 Ortholog",
         color = "Cluster Pair") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
    scale_y_continuous(limits = c(0, 1)) + # Set x-axis limits
    geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
    annotate("text", x = 0.2, y = 0.95, label = "0.2", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
    annotate("text", x = 0.5, y = 0.95, label = "0.5", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
    annotate("text", x = 1, y = 0.95, label = "1", hjust = -0.2, vjust = -0.5) +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}


#' PlotSubclassDEIntersectionCDF
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the de-curves family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param all_genes1 (auto) parameter
#' @param all_genes2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param output_path (auto) parameter
#' @param cluster_pairs (auto) parameter
#' @param pair_colors (auto) parameter
#' @param normalize.within (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-curves
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassDEIntersectionCDF <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, output_path, cluster_pairs, pair_colors, normalize.within = TRUE) {
  # Extract dataframes
  df1 <- list1$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  df2 <- list2$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  # Calculate cumulative distribution for each specified cluster pair
  for (i in seq_along(cluster_pairs)) {
    pair <- cluster_pairs[[i]]
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1) %>% filter(gene %in% intersecting_genes)
    df2_filtered <- df2 %>% filter(cluster == cluster2) %>% filter(gene %in% intersecting_genes)
    
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(c(df1_filtered$avg_log2FC, df2_filtered$avg_log2FC), na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 50)
    
    total_intersecting_genes_init <- length(union(df1_filtered$gene, df2_filtered$gene))
    
    count <- c()
    for (l in log2FC_grid) {
      de_genes_1 <- df1_filtered$gene[df1_filtered$avg_log2FC > l]
      de_genes_2 <- df2_filtered$gene[df2_filtered$avg_log2FC > l]
      intersecting_de_genes <- intersect(de_genes_1, de_genes_2)
      total_intersecting_genes <- length(union(de_genes_1, de_genes_2))
      if (normalize.within == TRUE) { 
        count <- c(count, length(intersecting_de_genes) / total_intersecting_genes) }
      else { count <- c(count, length(intersecting_de_genes) / total_intersecting_genes_init) }
      
    }
    
    # Calculate the cumulative number of shared genes as avg_log2FC varies
    shared_gene_counts <<- data.frame(
      avg_log2FC = log2FC_grid,
      count = count,
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      Color = pair_colors[i]
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
    # # Save intersecting genes to a text file
    # filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_intersecting_genes.txt")
    # write.table(intersecting_de_genes, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  cdf_data$Cluster1 <- factor(cdf_data$Cluster1, levels = subclass.order[subclass.order %in% cdf_data$Cluster1])
  cdf_data$Cluster2 <- factor(cdf_data$Cluster2, levels = subclass.order[subclass.order %in% cdf_data$Cluster2])
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count, color = interaction(Cluster1, Cluster2), group = interaction(Cluster1, Cluster2))) +
    geom_line(size = 1) + # Use line size
    scale_color_manual(values = setNames(pair_colors, unique(interaction(cdf_data$Cluster1, cdf_data$Cluster2))), guide = guide_legend(reverse = TRUE)) + # Set the colors using pair_colors
    labs(title = paste0(""),
         x = "avg_log2FC",
         y = "Fraction of 1:1 DE Genes Shared Across Species",
         color = "Cluster Pair") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
    geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
    annotate("text", x = 0.2, y = 0.5, label = "0.2", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
    annotate("text", x = 0.5, y = 0.5, label = "0.5", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
    annotate("text", x = 1, y = 0.5, label = "1", hjust = -0.2, vjust = -0.5) +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}


#' PlotSubclassDEIntersectionCDFTopX
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the de-curves family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param all_genes1 (auto) parameter
#' @param all_genes2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param output_path (auto) parameter
#' @param cluster_pairs (auto) parameter
#' @param pair_colors (auto) parameter
#' @param top_genes_seq (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-curves
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassDEIntersectionCDFTopX <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, output_path, cluster_pairs, pair_colors, top_genes_seq) {
  # Extract dataframes
  df1 <- list1$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  df2 <- list2$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  # Calculate cumulative distribution for each specified cluster pair
  for (i in seq_along(cluster_pairs)) {
    pair <- cluster_pairs[[i]]
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1) %>% filter(gene %in% intersecting_genes)
    df2_filtered <- df2 %>% filter(cluster == cluster2) %>% filter(gene %in% intersecting_genes)
    
    # Rank genes by avg_log2FC
    df1_filtered <- df1_filtered %>% arrange(dplyr::desc(avg_log2FC))
    df2_filtered <- df2_filtered %>% arrange(dplyr::desc(avg_log2FC))
    
    # Calculate the cumulative fraction of shared genes as top X DE genes vary
    shared_gene_counts <- data.frame(
      top_genes = top_genes_seq,
      count = sapply(top_genes_seq, function(x) {
        top_genes1 <- head(df1_filtered$gene, x)
        top_genes2 <- head(df2_filtered$gene, x)
        length(intersect(top_genes1, top_genes2)) / length(union(top_genes1, top_genes2))
      }),
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      Color = pair_colors[i]
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
  }
  
  cdf_data$Cluster1 <- factor(cdf_data$Cluster1, levels = subclass.order[subclass.order %in% cdf_data$Cluster1])
  cdf_data$Cluster2 <- factor(cdf_data$Cluster2, levels = subclass.order[subclass.order %in% cdf_data$Cluster2])
  cdf_data$ClusterPair <- interaction(cdf_data$Cluster1, cdf_data$Cluster2)
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = top_genes, y = count, color = ClusterPair, group = ClusterPair)) +
    geom_line(size = 1) + # Use line size
    scale_color_manual(values = setNames(pair_colors, unique(cdf_data$ClusterPair)), guide = guide_legend(reverse = TRUE)) + # Set the colors using pair_colors
    labs(title = paste0(""),
         x = "Top X DE Genes",
         y = "Fraction of 1:1 DE Genes Shared Across Species",
         color = "Cluster Pair") +
    scale_x_reverse() + # Reverse the x-axis
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}


#' PlotSubclassDEIntersectionScatter
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the de-curves family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param cluster_pairs (auto) parameter
#' @param pair_colors (auto) parameter
#' @param log2FC_threshold (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-curves
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassDEIntersectionScatter <- function(list1, list2, sample.name1, sample.name2, subclass.order, cluster_pairs, pair_colors, log2FC_threshold = 0.5) {
  # Extract dataframes
  df1 <- list1$subclass %>%
           filter(pct.1 >= 0.2) %>% 
           filter(p_val_adj < 0.05)
  df2 <- list2$subclass %>%
           filter(pct.1 >= 0.2) %>% 
           filter(p_val_adj < 0.05)
  
  # Initialize an empty data frame to store counts
  scatter_data <- data.frame()
  
  # Calculate DE gene counts for each specified cluster pair
  for (i in seq_along(cluster_pairs)) {
    pair <- cluster_pairs[[i]]
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Count DE genes above the specified log2FC threshold
    de_genes_1_count <- sum(df1_filtered$avg_log2FC > log2FC_threshold, na.rm = TRUE)
    de_genes_2_count <- sum(df2_filtered$avg_log2FC > log2FC_threshold, na.rm = TRUE)
    
    # Store the counts and cluster pair information
    scatter_data <- rbind(scatter_data, data.frame(
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      DE_Genes_List1 = de_genes_1_count,
      DE_Genes_List2 = de_genes_2_count,
      Color = pair_colors[i]
    ))
  }
  
  scatter_data$Cluster1 <- factor(scatter_data$Cluster1, levels = subclass.order[subclass.order %in% scatter_data$Cluster1])
  scatter_data$Cluster2 <- factor(scatter_data$Cluster2, levels = subclass.order[subclass.order %in% scatter_data$Cluster2])
  
  # Determine the range for the axes
  max_count <- max(c(scatter_data$DE_Genes_List1, scatter_data$DE_Genes_List2), na.rm = TRUE)
  
  # Plot the scatter plot with equal and square axes, and a diagonal line
  ggplot(scatter_data, aes(x = DE_Genes_List1, y = DE_Genes_List2, color = interaction(Cluster1, Cluster2))) +
    geom_point(size = 3) + # Adjust point size as needed
    scale_color_manual(values = setNames(pair_colors, unique(interaction(scatter_data$Cluster1, scatter_data$Cluster2))), guide = guide_legend(reverse = TRUE)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") + # Add diagonal line
    coord_fixed(ratio = 1) + # Ensure axes are equal and square
    scale_x_continuous(limits = c(0, max_count)) + 
    scale_y_continuous(limits = c(0, max_count)) +
    labs(title = paste0("log2FC > ", log2FC_threshold),
         x = paste0("Number of DE Genes in ", sample.name1),
         y = paste0("Number of DE Genes in ", sample.name2),
         color = "Cluster Pair") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}


#' PlotSubclassDEIntersectionOverlapScatter
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the de-curves family.
#' @param list1 (auto) parameter
#' @param list2 (auto) parameter
#' @param all_genes1 (auto) parameter
#' @param all_genes2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param cluster_pairs (auto) parameter
#' @param pair_colors (auto) parameter
#' @param log2FC_threshold (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-curves
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassDEIntersectionOverlapScatter <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, cluster_pairs, pair_colors, log2FC_threshold = 0.5) {
  # Extract dataframes
  df1 <- list1$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  df2 <- list2$subclass %>%
    filter(pct.1 >= 0.2) %>% 
    filter(p_val_adj < 0.05)
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Initialize an empty data frame to store counts
  scatter_data_1 <- data.frame()
  scatter_data_2 <- data.frame()
  
  # Calculate DE gene counts for each specified cluster pair
  for (i in seq_along(cluster_pairs)) {
    pair <- cluster_pairs[[i]]
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Count DE genes above the specified log2FC threshold
    de_genes_1 <- df1_filtered %>% filter(avg_log2FC > log2FC_threshold) %>% pull(gene)
    de_genes_2 <- df2_filtered %>% filter(avg_log2FC > log2FC_threshold) %>% pull(gene)
    
    de_genes_1_count <- sum(de_genes_1 %in% intersecting_genes)
    de_genes_2_count <- sum(de_genes_2 %in% intersecting_genes)
    
    # Count intersecting DE genes between the two lists
    intersecting_de_genes_count <- length(intersect(de_genes_1, de_genes_2))
    
    # Store the counts and cluster pair information for both lists
    scatter_data_1 <- rbind(scatter_data_1, data.frame(
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      DE_Genes_List = de_genes_1_count,
      Intersecting_DE_Genes = intersecting_de_genes_count,
      Color = pair_colors[i],
      Sample = sample.name1
    ))
    
    scatter_data_2 <- rbind(scatter_data_2, data.frame(
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      DE_Genes_List = de_genes_2_count,
      Intersecting_DE_Genes = intersecting_de_genes_count,
      Color = pair_colors[i],
      Sample = sample.name2
    ))
  }
  
  scatter_data_1$Cluster1 <- factor(scatter_data_1$Cluster1, levels = subclass.order[subclass.order %in% scatter_data_1$Cluster1])
  scatter_data_1$Cluster2 <- factor(scatter_data_1$Cluster2, levels = subclass.order[subclass.order %in% scatter_data_1$Cluster2])
  
  scatter_data_2$Cluster1 <- factor(scatter_data_2$Cluster1, levels = subclass.order[subclass.order %in% scatter_data_2$Cluster1])
  scatter_data_2$Cluster2 <- factor(scatter_data_2$Cluster2, levels = subclass.order[subclass.order %in% scatter_data_2$Cluster2])
  
  # Determine the maximum range for the axes
  max_x <- max(c(scatter_data_1$DE_Genes_List, scatter_data_2$DE_Genes_List), na.rm = TRUE)
  max_y <- max(c(scatter_data_1$Intersecting_DE_Genes, scatter_data_2$Intersecting_DE_Genes), na.rm = TRUE)
  max_limit <- max(max_x, max_y)
  
  # Plot the scatter plots for both lists
  plot1 <- ggplot(scatter_data_1, aes(x = DE_Genes_List, y = Intersecting_DE_Genes, color = interaction(Cluster1, Cluster2))) +
    geom_point(size = 3) + # Adjust point size as needed
    scale_color_manual(values = setNames(pair_colors, unique(interaction(scatter_data_1$Cluster1, scatter_data_1$Cluster2))), guide = guide_legend(reverse = TRUE)) +
    coord_fixed(ratio = 1) + # Ensure axes are equal and square
    scale_x_continuous(limits = c(0, max_limit)) + 
    scale_y_continuous(limits = c(0, max_limit)) +
    labs(title = paste0("log2FC > ", log2FC_threshold, ": ", sample.name1),
         x = paste0("Number of DE Genes in ", sample.name1),
         y = "Number of Intersecting DE Genes",
         color = "Cluster Pair") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  plot2 <- ggplot(scatter_data_2, aes(x = DE_Genes_List, y = Intersecting_DE_Genes, color = interaction(Cluster1, Cluster2))) +
    geom_point(size = 3) + # Adjust point size as needed
    scale_color_manual(values = setNames(pair_colors, unique(interaction(scatter_data_2$Cluster1, scatter_data_2$Cluster2))), guide = guide_legend(reverse = TRUE)) +
    coord_fixed(ratio = 1) + # Ensure axes are equal and square
    scale_x_continuous(limits = c(0, max_limit)) + 
    scale_y_continuous(limits = c(0, max_limit)) +
    labs(title = paste0("log2FC > ", log2FC_threshold, ": ", sample.name2),
         x = paste0("Number of DE Genes in ", sample.name2),
         y = "Number of Intersecting DE Genes",
         color = "Cluster Pair") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  # Return the plots as a list
  list(plot1 = plot1, plot2 = plot2)
}


#' PlotSubclassGeneCountCDFDiff
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the de-curves family.
#' @param de_df_1 (auto) parameter
#' @param de_df_2 (auto) parameter
#' @param subclass_pairs (auto) parameter
#' @param subclass_colors (auto) parameter
#' @param min.pct (auto) parameter
#' @param max.pval (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-curves
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassGeneCountCDFDiff <- function(de_df_1, de_df_2, subclass_pairs, subclass_colors, min.pct = 0.1, max.pval = 0.05) {
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  for (pair in subclass_pairs) {
    sbcl_1 <- pair[1]
    sbcl_2 <- pair[2]
    
    # Filter dataframes for the specific subclass
    df_filtered_1 <- de_df_1$subclass %>% filter(cluster == sbcl_1, pct.1 >= min.pct, p_val_adj < max.pval)
    df_filtered_2 <- de_df_2$subclass %>% filter(cluster == sbcl_2, pct.1 >= min.pct, p_val_adj < max.pval)
    
    # Calculate the max avg_log2FC for each group
    max_log2FC_1 <- max(df_filtered_1$avg_log2FC, na.rm = TRUE)
    max_log2FC_2 <- max(df_filtered_2$avg_log2FC, na.rm = TRUE)
    
    # Create a common log2FC grid
    log2FC_grid <- seq(0.2, min(max(max_log2FC_1, max_log2FC_2), 2), length.out = 50)
    
    count_1 <- sapply(log2FC_grid, function(l) sum(df_filtered_1$avg_log2FC > l, na.rm = TRUE))
    count_2 <- sapply(log2FC_grid, function(l) sum(df_filtered_2$avg_log2FC > l, na.rm = TRUE))
    
    # Calculate the difference in the number of DE genes
    diff_counts <- count_1 - count_2
    
    # Create data frame
    diff_gene_counts <- data.frame(
      avg_log2FC = log2FC_grid,
      count_difference = diff_counts,
      Comparison = paste(sbcl_1, "vs", sbcl_2)
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, diff_gene_counts)
  }
  
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count_difference, color = Comparison, group = Comparison)) +
    geom_line(size = 1) +
    scale_color_manual(values = subclass_colors) +
    labs(title = paste0("Comparison of DE Gene Counts"),
         x = "avg_log2FC",
         y = "Difference in Number of DE Genes",
         color = "Subclass Comparison") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.2, 2)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}


#' PlotPairwiseSubclassGeneCountCDFDifference
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the de-curves family.
#' @param de_results_list_1 (auto) parameter
#' @param de_results_list_2 (auto) parameter
#' @param subclass_pairs (auto) parameter
#' @param subclass_colors (auto) parameter
#' @param min.pct (auto) parameter
#' @param max.pval (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family de-curves
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotPairwiseSubclassGeneCountCDFDifference <- function(de_results_list_1, de_results_list_2, subclass_pairs, subclass_colors, min.pct = 0.1, max.pval = 0.05) {
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  for (pair in subclass_pairs) {
    sbcl_1_1 <- pair[[1]][1]
    sbcl_1_2 <- pair[[1]][2]
    sbcl_2_1 <- pair[[2]][1]
    sbcl_2_2 <- pair[[2]][2]
    
    de_results_1 <- de_results_list_1[[paste(sbcl_1_1, sbcl_1_2, sep = "_vs_")]]
    de_results_2 <- de_results_list_2[[paste(sbcl_2_1, sbcl_2_2, sep = "_vs_")]]
    
    de_results_1 <- de_results_1[de_results_1$pct.1 > min.pct & de_results_1$p_val_adj < max.pval, ]
    de_results_2 <- de_results_2[de_results_2$pct.1 > min.pct & de_results_2$p_val_adj < max.pval, ]
    
    # Separate DE genes for each subclass
    df_filtered_1 <<- de_results_1 # [de_results_1$avg_log2FC > 0, ]
    df_filtered_2 <<- de_results_2 # [de_results_2$avg_log2FC > 0, ]
    
    # Calculate the max avg_log2FC for each group
    max_log2FC_1 <- max(df_filtered_1$avg_log2FC, na.rm = TRUE)
    max_log2FC_2 <- max(df_filtered_2$avg_log2FC, na.rm = TRUE)
    
    # Create a common log2FC grid
    log2FC_grid <- seq(0.2, min(max(max_log2FC_1, max_log2FC_2), 2), length.out = 50)
    
    count_1 <- sapply(log2FC_grid, function(l) sum(df_filtered_1$avg_log2FC > l, na.rm = TRUE))
    count_2 <- sapply(log2FC_grid, function(l) sum(df_filtered_2$avg_log2FC > l, na.rm = TRUE))
    
    # Calculate the difference in the number of DE genes between species
    diff_counts <- count_1 - count_2
    
    # Create data frame
    diff_gene_counts <- data.frame(
      avg_log2FC = log2FC_grid,
      count_difference = diff_counts,
      Comparison = paste(sbcl_1_1, "vs", sbcl_1_2, "minus", sbcl_2_1, "vs", sbcl_2_2)
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, diff_gene_counts)
  }
  
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count_difference, color = Comparison, group = Comparison)) +
    geom_line(size = 1) +
    scale_color_manual(values = subclass_colors) +
    labs(title = paste0("Pairwise Comparison of DE Gene Counts between Species"),
         x = "avg_log2FC",
         y = "Difference in Number of DE Genes",
         color = "Subclass Comparison") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.49, 2)) +
    scale_y_continuous(limits = c(-150, 150)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}
