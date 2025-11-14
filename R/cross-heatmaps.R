#' IntegratedClusterOverlapHeatmap
#'
#' Generate overlap heatmap showing shared cluster membership between integrated datasets.
#' 
#' This function quantifies how cells from different datasets (e.g., different species,
#' conditions, or batches) co-cluster in integrated space. For each pair of cell type
#' identities across datasets, it calculates the overlap coefficient: the minimum of the
#' two fractions of cells sharing the same integrated cluster. This metric is robust to
#' differences in cell type abundance between datasets and emphasizes strong co-clustering
#' relationships.
#' 
#' @param integrated.obj Integrated Seurat object containing cells from multiple datasets.
#'   Must have metadata columns specified by `integvar.col`, `ident.col`, and `cluster.col`.
#' @param integvar.col String specifying the metadata column that distinguishes datasets
#'   (e.g., "species", "condition", "batch"). Should have exactly two unique values.
#' @param ident.col String specifying the metadata column containing cell type identities
#'   to compare (e.g., "subclass", "type"). These labels should be assigned independently
#'   in each dataset prior to integration.
#' @param cluster.col String specifying the metadata column containing integrated cluster
#'   assignments (e.g., "integrated_snn_res.1"). These clusters are computed from the
#'   integrated space.
#' @param primary_order_row Character vector specifying the desired order of cell type
#'   labels for the first dataset (y-axis). Used to match and sort identities; see
#'   \code{\link{sort_idents}} for pattern matching behavior.
#' @param primary_order_col Character vector specifying the desired order of cell type
#'   labels for the second dataset (x-axis). Typically identical to `primary_order_row`
#'   unless datasets use different naming schemes.
#' @param col.low String specifying the color for low overlap values. Default is "white".
#' @param col.high String specifying the color for high overlap values. Default is "red".
#' @param x.lab.rot Logical indicating whether to rotate x-axis labels 90 degrees. Default
#'   is TRUE, useful for long cell type names.
#' @param show_text Logical indicating whether to display percentage values in each tile.
#'   Default is TRUE. Set to FALSE for cleaner visualizations when categories are numerous.
#'   
#' @details
#' **Overlap calculation**:
#' 
#' For each pair of cell types (A from dataset 1, B from dataset 2), the function:
#' 
#' 1. Iterates over all integrated clusters
#' 2. For each cluster k, calculates:
#'    - f_A,k = fraction of cells in type A that belong to cluster k
#'    - f_B,k = fraction of cells in type B that belong to cluster k
#' 3. Computes the minimum: min(f_A,k, f_B,k)
#' 4. Sums these minima across all clusters to get the total overlap
#' 
#' This overlap coefficient ranges from 0 (no co-clustering) to 1 (perfect co-clustering).
#' Values are displayed as percentages (0-100%).
#' 
#' **Row ordering**:
#' Cell types are sorted using the \code{\link{sort_idents}} function, which handles
#' complex naming patterns (e.g., "IT_A_1", "IT_A_10") and applies the user-specified
#' primary order. Rows are reversed in the final plot (last identity at top) to match
#' common dendrogram orientations.
#' 
#' @return A ggplot2 heatmap object showing overlap percentages between cell types across
#'   datasets. The plot includes axis labels indicating dataset names and uses a white-to-red
#'   gradient (or custom colors) to represent overlap strength.
#' 
#' @export
#' @family integration-heatmaps
#' 
#' @examples
#' \dontrun{
#'  # Basic usage with subclass-level comparison
#'  subclass_order <- c("IT", "L5PT", "L6CT", "L6b", "Pvalb", "Sst", 
#'                      "Vip", "Lamp5", "Astro", "Oligo", "OPC")
#'  
#'  p <- IntegratedClusterOverlapHeatmap(
#'    obj_integrated,
#'    integvar.col = "species",
#'    ident.col = "subclass",
#'    cluster.col = "integrated_snn_res.1",
#'    primary_order_row = subclass_order,
#'    primary_order_col = subclass_order
#'  )
#'  print(p)
#'  
#'  # Different orders for species with different nomenclature
#'  IntegratedClusterOverlapHeatmap(
#'    obj_integrated,
#'    integvar.col = "species",
#'    ident.col = "subclass",
#'    cluster.col = "integrated_snn_res.1",
#'    primary_order_row = c("L2/3", "L4", "L5IT", "L6IT"),  # Mouse names
#'    primary_order_col = c("IT_A", "IT_B", "IT_C", "IT_D")  # Opossum names
#'  )
#'  
#'  # Cleaner visualization without text labels
#'  IntegratedClusterOverlapHeatmap(
#'    obj_integrated,
#'    integvar.col = "condition",
#'    ident.col = "type",
#'    cluster.col = "integrated_snn_res.2",
#'    primary_order_row = type_order,
#'    primary_order_col = type_order,
#'    show_text = FALSE,
#'    col.high = "darkblue"
#'  )
#' }
IntegratedClusterOverlapHeatmap <- function(integrated.obj, integvar.col, ident.col, cluster.col, 
                                            primary_order_row, primary_order_col, 
                                            col.low = "white", col.high = "red", 
                                            x.lab.rot = TRUE, show_text = TRUE) {
  
  # Extract the relevant columns from metadata
  metadata <- integrated.obj@meta.data %>% select(all_of(c(integvar.col, cluster.col, ident.col)))
  
  # Rename columns for convenience
  colnames(metadata) <- c("integvar", "cluster", "ident")
  
  # Get unique integvar levels (should be exactly 2 for pairwise comparison)
  integvar_levels <- unique(metadata$integvar)
  
  # Initialize an empty list to store overlap matrices for each pair of integvar levels
  overlap_matrices <- list()
  
  # Sorting function with separate primary orders for rows and columns
  sort_ident <- function(ident, primary_order) {
    # Sort primary_order by decreasing length to match longer patterns first
    primary_order_sorted <- primary_order[order(nchar(primary_order), decreasing = TRUE)]
    
    # Extract primary identifier by matching against primary_order
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order_sorted, collapse = "|")))
    
    # Extract suffix after the last underscore
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    
    # Convert numeric suffixes for proper numerical sorting
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    
    # For non-numeric suffixes, prefix with "Z" to sort after numeric suffixes
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    
    # Create sorting dataframe
    df <- data.frame(
      ident = ident, 
      primary = primary, 
      suffix = suffix, 
      suffix_numeric = suffix_numeric,
      stringsAsFactors = FALSE
    )
    
    # Sort by: (1) position in primary_order, (2) numeric suffix, (3) alphabetic suffix
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    
    return(df$ident)
  }
  
  # Loop through each pair of integvar levels
  for (i in 1:(length(integvar_levels) - 1)) {
    for (j in (i + 1):length(integvar_levels)) {
      integvar1 <- integvar_levels[i]
      integvar2 <- integvar_levels[j]
      
      # Filter the metadata for the two integvar levels
      data_integvar1 <- metadata %>% filter(integvar == integvar1)
      data_integvar2 <- metadata %>% filter(integvar == integvar2)
      
      # Get unique ident levels for each integvar level
      ident_levels1 <- unique(data_integvar1$ident)
      ident_levels2 <- unique(data_integvar2$ident)
      
      # Sort ident levels separately for rows and columns
      sorted_ident_levels1 <- sort_ident(ident_levels1, primary_order_row)
      sorted_ident_levels2 <- sort_ident(ident_levels2, primary_order_col)
      
      # Get unique clusters
      clusters <- unique(metadata$cluster)
      
      # Initialize the overlap matrix
      # Note: rows are reversed for display (last identity at top)
      overlap_matrix <- matrix(
        0, 
        nrow = length(sorted_ident_levels1), 
        ncol = length(sorted_ident_levels2), 
        dimnames = list(rev(sorted_ident_levels1), sorted_ident_levels2)
      )
      
      # Calculate overlap fractions using minimum overlap coefficient
      for (ident1 in sorted_ident_levels1) {
        for (ident2 in sorted_ident_levels2) {
          for (cluster in clusters) {
            # Fraction of ident1 cells in this cluster
            n_ident1_in_cluster <- sum(data_integvar1$ident == ident1 & data_integvar1$cluster == cluster)
            n_ident1_total <- sum(data_integvar1$ident == ident1)
            fraction_ident1 <- if (n_ident1_total > 0) n_ident1_in_cluster / n_ident1_total else 0
            
            # Fraction of ident2 cells in this cluster
            n_ident2_in_cluster <- sum(data_integvar2$ident == ident2 & data_integvar2$cluster == cluster)
            n_ident2_total <- sum(data_integvar2$ident == ident2)
            fraction_ident2 <- if (n_ident2_total > 0) n_ident2_in_cluster / n_ident2_total else 0
            
            # Add minimum of the two fractions (overlap coefficient)
            overlap_matrix[ident1, ident2] <- overlap_matrix[ident1, ident2] + min(fraction_ident1, fraction_ident2)
          }
        }
      }
      
      # Convert overlap fractions to percentages
      overlap_matrix <- overlap_matrix * 100
      
      # Store the overlap matrix in the list
      overlap_matrices[[paste(integvar1, integvar2, sep = "_vs_")]] <- overlap_matrix
    }
  }
  
  # Plot the heatmap for each pair of integvar levels
  # Note: typically only one pair for two-dataset comparisons
  for (name in names(overlap_matrices)) {
    overlap_matrix <- overlap_matrices[[name]]
    
    # Melt the matrix for ggplot
    melted <- melt(overlap_matrix)
    colnames(melted) <- c("row", "col", "Percentage")
    
    # Create the heatmap plot
    p <- ggplot(melted, aes(y = row, x = col)) + 
      geom_tile(aes(fill = Percentage)) + 
      scale_fill_gradient(low = col.low, high = col.high, limits = c(0, 100)) + 
      theme_bw() + 
      xlab(integvar_levels[2]) + 
      ylab(integvar_levels[1]) + 
      theme(
        axis.text.x = element_text(size = 16, face = "italic", hjust = 1, angle = ifelse(x.lab.rot, 90, 0)),
        axis.text.y = element_text(size = 16, face = "italic"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
      ) +
      coord_fixed()
    
    # Add text labels if show_text is TRUE
    if (show_text) {
      p <- p + geom_text(aes(label = sprintf("%.0f", Percentage)), size = 5)
    }
    
    # Ensure correct rotation for x-axis labels
    if (x.lab.rot) {
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    } else {
      p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
    }
    
    # Add title
    p <- p + ggtitle(paste(name, ident.col, "at", cluster.col))
  }
  
  return(p)
}


#' IntegratedClusterMakeupHeatmap
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the integration-heatmaps family.
#' @param integrated.obj (auto) parameter
#' @param integvar.col (auto) parameter
#' @param ident.col (auto) parameter
#' @param cluster.col (auto) parameter
#' @param primary_order (auto) parameter
#' @param col.low (auto) parameter
#' @param col.high (auto) parameter
#' @param x.lab.rot (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family integration-heatmaps
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
IntegratedClusterMakeupHeatmap <- function(integrated.obj, integvar.col, ident.col, cluster.col, primary_order, 
                                           col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  # Extract the relevant columns
  metadata <- integrated.obj@meta.data %>% select(all_of(c(integvar.col, cluster.col, ident.col)))
  
  # Rename columns for convenience
  colnames(metadata) <- c("integvar", "cluster", "ident")
  
  # Get unique integvar levels
  integvar_levels <- unique(metadata$integvar)
  
  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(df$ident)
  }
  
  # Initialize a list to store plots
  plot_list <- list()
  cluster_levels <- unique(metadata$cluster)
  cluster_levels <- sort(cluster_levels)
  
  # Plot the fraction of each cluster.col that comes from each ident.col, split by integvar
  for (iv in integvar_levels) {
    data_integvar <- metadata %>% filter(integvar == iv)
    ident_levels <- unique(data_integvar$ident)
    sorted_ident_levels <- sort_ident(ident_levels, primary_order)

    fraction_matrix <- matrix(0, nrow = length(cluster_levels), ncol = length(sorted_ident_levels), 
                              dimnames = list(cluster_levels, sorted_ident_levels))
    
    for (cluster in cluster_levels) {
      for (ident in sorted_ident_levels) {
        if (any(data_integvar$cluster == cluster)) {
        fraction_matrix[cluster, ident] <- sum(data_integvar$cluster == cluster & data_integvar$ident == ident) / sum(metadata$cluster == cluster)
        # fraction_matrix[cluster, ident] <- sum(data_integvar$cluster == cluster & data_integvar$ident == ident) / sum(data_integvar$ident == ident)
        }
        else { fraction_matrix[cluster, ident] <- 0 }
      }
    }
    
    # Melt the matrix for ggplot
    fraction_matrix <- fraction_matrix * 100
    melted_fraction <- melt(fraction_matrix)
    colnames(melted_fraction) <- c("row", "col", "Fraction")
    
    # Create the heatmap plot
    p_fraction <- ggplot(melted_fraction, aes(x = col, y = factor(row, levels = rev(cluster_levels)))) + 
      geom_tile(aes(fill = Fraction)) + 
      scale_fill_gradient(low = col.low, high = col.high, limits = c(0, 100)) + 
      geom_text(aes(label = sprintf("%.0f", Fraction)), size = 3) +
      theme_bw() + 
      xlab(iv) + 
      ylab("") + 
      theme(axis.text.x = element_text(size=12, face="italic", hjust=1, angle = ifelse(x.lab.rot, 90, 0)), 
            axis.text.y = element_text(size=12, face="italic"), 
            axis.title.x = element_text(size=12), 
            axis.title.y = element_text(size=12)) +
      coord_fixed(ratio = length(ident_levels) / length(cluster_levels))
    
    # Ensure correct rotation
    if (x.lab.rot) {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    } else {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    }
    
    plot_list[[iv]] <- p_fraction + theme(legend.position = "none")
    if (iv == integvar_levels[1]) { plot_list[[iv]] <- plot_list[[iv]] + ggtitle("% Makeup of Integrated Clusters") }
    else { plot_list[[iv]] <- plot_list[[iv]] + ggtitle("") }
  }
  return(plot_list)
}


#' IdentToIntegratedClusterHeatmap
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the integration-heatmaps family.
#' @param integrated.obj (auto) parameter
#' @param integvar.col (auto) parameter
#' @param ident.col (auto) parameter
#' @param cluster.col (auto) parameter
#' @param primary_order (auto) parameter
#' @param col.low (auto) parameter
#' @param col.high (auto) parameter
#' @param x.lab.rot (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family integration-heatmaps
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
IdentToIntegratedClusterHeatmap <- function(integrated.obj, integvar.col, ident.col, cluster.col, primary_order, 
                                           col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  # Extract the relevant columns
  metadata <- integrated.obj@meta.data %>% select(all_of(c(integvar.col, cluster.col, ident.col)))
  
  # Rename columns for convenience
  colnames(metadata) <- c("integvar", "cluster", "ident")
  
  # Get unique integvar levels
  integvar_levels <- unique(metadata$integvar)
  
  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(df$ident)
  }
  
  # Initialize a list to store plots
  plot_list <- list()
  cluster_levels <- unique(metadata$cluster)
  cluster_levels <- sort(cluster_levels)
  
  # Plot the fraction of each cluster.col that comes from each ident.col, split by integvar
  for (iv in integvar_levels) {
    data_integvar <- metadata %>% filter(integvar == iv)
    ident_levels <- unique(data_integvar$ident)
    sorted_ident_levels <- sort_ident(ident_levels, primary_order)
    
    fraction_matrix <- matrix(0, nrow = length(cluster_levels), ncol = length(sorted_ident_levels), 
                              dimnames = list(cluster_levels, sorted_ident_levels))
    
    for (cluster in cluster_levels) {
      for (ident in sorted_ident_levels) {
        if (any(data_integvar$cluster == cluster)) {
          fraction_matrix[cluster, ident] <- sum(data_integvar$cluster == cluster & data_integvar$ident == ident) / sum(data_integvar$ident == ident)
        }
        else { fraction_matrix[cluster, ident] <- 0 }
      }
    }
    
    # Melt the matrix for ggplot
    fraction_matrix <- fraction_matrix * 100
    melted_fraction <- melt(fraction_matrix)
    colnames(melted_fraction) <- c("row", "col", "Fraction")
    
    # Create the heatmap plot
    p_fraction <- ggplot(melted_fraction, aes(x = col, y = factor(row, levels = rev(cluster_levels)))) + 
      geom_tile(aes(fill = Fraction)) + 
      scale_fill_gradient(low = col.low, high = col.high, limits = c(0, 100)) + 
      geom_text(aes(label = sprintf("%.0f", Fraction)), size = 3) +
      theme_bw() + 
      xlab(iv) + 
      ylab("") + 
      theme(axis.text.x = element_text(size=12, face="italic", hjust=1, angle = ifelse(x.lab.rot, 90, 0)), 
            axis.text.y = element_text(size=12, face="italic"), 
            axis.title.x = element_text(size=12), 
            axis.title.y = element_text(size=12)) +
      coord_fixed(ratio = length(ident_levels) / length(cluster_levels))
    
    # Ensure correct rotation
    if (x.lab.rot) {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    } else {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    }
    
    plot_list[[iv]] <- p_fraction + theme(legend.position = "none")
    if (iv == integvar_levels[1]) { plot_list[[iv]] <- plot_list[[iv]] + ggtitle(paste("% of", ident.col, "in Integrated Clusters")) }
    else { plot_list[[iv]] <- plot_list[[iv]] + ggtitle("") }
  }
  return(plot_list)
}
