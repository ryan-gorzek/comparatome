#' IntegratedClusterOverlapHeatmap
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the integration-heatmaps family.
#' @param integrated.obj (auto) parameter
#' @param integvar.col (auto) parameter
#' @param ident.col (auto) parameter
#' @param cluster.col (auto) parameter
#' @param primary_order_row (auto) parameter
#' @param primary_order_col (auto) parameter
#' @param col.low (auto) parameter
#' @param col.high (auto) parameter
#' @param x.lab.rot (auto) parameter
#' @param show_text (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family integration-heatmaps
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
IntegratedClusterOverlapHeatmap <- function(integrated.obj, integvar.col, ident.col, cluster.col, 
                                            primary_order_row, primary_order_col, 
                                            col.low = "white", col.high = "red", 
                                            x.lab.rot = TRUE, show_text = TRUE) {
  
  # Extract the relevant columns
  metadata <- integrated.obj@meta.data %>% select(all_of(c(integvar.col, cluster.col, ident.col)))
  
  # Rename columns for convenience
  colnames(metadata) <- c("integvar", "cluster", "ident")
  
  # Get unique integvar levels
  integvar_levels <- unique(metadata$integvar)
  
  # Initialize an empty list to store overlap matrices for each pair of integvar levels
  overlap_matrices <- list()
  
  # Sorting function with separate primary orders for rows and columns
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
      overlap_matrix <- matrix(0, nrow = length(sorted_ident_levels1), ncol = length(sorted_ident_levels2), 
                               dimnames = list(rev(sorted_ident_levels1), sorted_ident_levels2))
      
      # Calculate overlap fractions
      for (ident1 in sorted_ident_levels1) {
        for (ident2 in sorted_ident_levels2) {
          for (cluster in clusters) {
            fraction_ident1 <- sum(data_integvar1$ident == ident1 & data_integvar1$cluster == cluster) / sum(data_integvar1$ident == ident1)
            fraction_ident2 <- sum(data_integvar2$ident == ident2 & data_integvar2$cluster == cluster) / sum(data_integvar2$ident == ident2)
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
  for (name in names(overlap_matrices)) {
    overlap_matrix <- overlap_matrices[[name]]
    
    # Melt the matrix for ggplot
    melted <- reshape2::melt(overlap_matrix)
    colnames(melted) <- c("row", "col", "Percentage")
    
    # Create the heatmap plot
    p <- ggplot(melted, aes(y = row, x = col)) + 
      geom_tile(aes(fill = Percentage)) + 
      scale_fill_gradient(low = col.low, high = col.high, limits = c(0, 100)) + 
      theme_bw() + 
      xlab(integvar_levels[2]) + 
      ylab(integvar_levels[1]) + 
      theme(axis.text.x = element_text(size=16, face="italic", hjust=1, angle = ifelse(x.lab.rot, 90, 0)),
            axis.text.y = element_text(size=16, face="italic"),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16)) +
      coord_fixed()
    
    # Add text labels if show_text is TRUE
    if (show_text) {
      p <- p + geom_text(aes(label = sprintf("%.0f", Percentage)), size = 5)
    }
    
    # Ensure correct rotation
    if (x.lab.rot) {
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    } else {
      p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
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
