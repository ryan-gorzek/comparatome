#' PlotPCLoadingsCorrelation
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the components family.
#' @param seurat_objects (auto) parameter
#' @param object_names (auto) parameter
#' @param num_pcs (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family components
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotPCLoadingsCorrelation <- function(seurat_objects, object_names, num_pcs = 10) {
  # Extract PCA loadings for the specified number of PCs
  loadings_list <- lapply(seurat_objects, function(obj) {
    Loadings(obj, reduction = "pca")[, 1:num_pcs]
  })
  
  # Find the intersection of genes between the two Seurat objects
  common_genes <- intersect(rownames(loadings_list[[1]]), rownames(loadings_list[[2]]))
  
  # Subset the loadings matrices to the common genes
  loadings_list <- lapply(loadings_list, function(loadings) {
    loadings[common_genes, ]
  })
  
  # Calculate the correlation matrix between the two sets of loadings
  cor_matrix <- cor(loadings_list[[1]], loadings_list[[2]])
  
  # Convert correlation matrix to long format for ggplot2
  cor_df <- melt(cor_matrix)
  colnames(cor_df) <- c("PC1", "PC2", "Correlation")
  
  # Create the heatmap using ggplot2
  ggplot(cor_df, aes(x = PC1, y = PC2, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = paste("PC Loadings Correlation of", object_names[1], "with", object_names[2]),
         x = paste0(object_names[1], " PC"),
         y = paste0(object_names[2], " PC")) +
    scale_y_discrete(limits = rev(levels(cor_df$PC1))) +
    coord_fixed() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' PlotVarianceExplainedWithinSpecies
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the components family.
#' @param seurat_obj (auto) parameter
#' @param species_names (auto) parameter
#' @param num_pcs (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family components
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotVarianceExplainedWithinSpecies <- function(seurat_obj, species_names, num_pcs = 10) {
  library(ggplot2)
  library(matrixStats)
  
  # Get PCA embeddings
  pca_embeddings <- Embeddings(seurat_obj, reduction = "pca")
  
  # Add PCA embeddings to metadata
  seurat_obj <- AddMetaData(seurat_obj, pca_embeddings, col.name = paste0("PC", 1:ncol(pca_embeddings)))
  
  # Initialize a dataframe to hold the variance explained data
  variance_data <- data.frame(PC = integer(), Variance_Explained = numeric(), Species = character())
  
  # Loop through each species to calculate the variance explained
  for (species_name in species_names) {
    # Subset the cells for the species
    cells <- WhichCells(seurat_obj, expression = species == species_name)
    
    # Fetch PCA data for the species
    pca_data <- FetchData(seurat_obj, vars = paste0("PC", 1:num_pcs), cells = cells)
    
    # Calculate the total variance using rowVars on the scale.data for this species
    species_data <- seurat_obj@assays$SCT@scale.data[, cells]
    total_variance <- sum(rowVars(species_data))
    
    # Calculate the variance for each PC within the species
    variances <- apply(pca_data, 2, sd)
    
    # Normalize by the total variance within the species
    variances_explained <- variances / total_variance * 100
    
    # Create a temporary dataframe for this species
    temp_df <- data.frame(
      PC = 1:num_pcs,
      Variance_Explained = variances_explained,
      Species = species_name
    )
    
    # Add this data to the main dataframe
    variance_data <- rbind(variance_data, temp_df)
  }

  # Plot the scree plot
  ggplot(variance_data, aes(x = PC, y = Variance_Explained, color = Species, group = Species)) +
    geom_line() +
    geom_point() +
    labs(title = paste("Within-Species Variance Explained by Top", num_pcs, "PCs"),
         x = "Principal Component",
         y = "% Variance Explained") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "red")) +
    scale_x_continuous(breaks = 1:num_pcs, labels = 1:num_pcs)
}


#' PlotCCAFeatureLoadings
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the components family.
#' @param seurat_obj1 (auto) parameter
#' @param seurat_obj2 (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family components
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotCCAFeatureLoadings <- function(seurat_obj1, seurat_obj2) {
  # Ensure both Seurat objects contain the same genes
  common_genes <- intersect(rownames(seurat_obj1), rownames(seurat_obj2))
  seurat_obj1 <- subset(seurat_obj1, features = common_genes)
  seurat_obj2 <- subset(seurat_obj2, features = common_genes)
  
  # Run CCA
  combined <- RunCCA(seurat_obj1, seurat_obj2, features = common_genes, rescale = TRUE)
  
  # Extract feature loadings for each gene in each object
  # The `Loadings` function retrieves feature loadings for the reduction
  feature_loadings <- Loadings(combined, reduction = "cca")
  
  # Convert the loadings to a data frame for easy plotting
  loadings_df <- as.data.frame(feature_loadings)
  loadings_df$Gene <- rownames(loadings_df)  # Add gene names for reference
  loadings_df <<- loadings_df
  # Plot the feature loadings (for the first two CCA components as an example)
  ggplot(loadings_df, aes(x = CC_1, y = CC_2)) +
    geom_point(alpha = 0.6) +
    labs(
      title = "CCA Feature Loadings Plot",
      x = "CC1",
      y = "CC2"
    ) +
    theme_minimal()
}


#' PlotElbowComparison
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the components family.
#' @param obj.mouse (auto) parameter
#' @param obj.opossum (auto) parameter
#' @param num.pcs (auto) parameter
#' @param colors (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family components
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotElbowComparison <- function(obj.mouse, obj.opossum, num.pcs = 30, colors = c("blue", "red")) {
  
  # Compute variance explained for mouse
  stdev.mouse <- obj.mouse@reductions$pca@stdev
  var.explained.mouse <- (stdev.mouse^2) / sum(stdev.mouse^2) * 100
  cum.var.mouse <- cumsum(var.explained.mouse)  # Cumulative VE
  df.mouse <- data.frame(PC = seq_along(var.explained.mouse), Variance = var.explained.mouse, 
                         CumulativeVariance = cum.var.mouse, Species = "Mouse")
  
  # Compute variance explained for opossum
  stdev.opossum <- obj.opossum@reductions$pca@stdev
  var.explained.opossum <- (stdev.opossum^2) / sum(stdev.opossum^2) * 100
  cum.var.opossum <- cumsum(var.explained.opossum)  # Cumulative VE
  df.opossum <- data.frame(PC = seq_along(var.explained.opossum), Variance = var.explained.opossum, 
                           CumulativeVariance = cum.var.opossum, Species = "Opossum")
  
  # Combine data
  df <- rbind(df.mouse, df.opossum)
  
  # Limit PCs for plotting
  df <- df[df$PC <= num.pcs, ]
  
  # Define a scaling factor to align the right y-axis (cumulative VE) with the left y-axis
  max_var <- max(df$Variance)  # Max single-PC variance explained
  max_cum_var <- max(df$CumulativeVariance)  # Max cumulative variance explained
  scale_factor <- max_var / max_cum_var  # Scale factor to align axes
  
  # Create the elbow plot with dual y-axes
  p <- ggplot(df, aes(x = PC, group = Species)) +
    # Left y-axis: Variance explained per PC (solid line)
    geom_line(aes(y = Variance, color = Species), size = 1) +  
    geom_point(aes(y = Variance, color = Species), size = 2) +
    
    # Right y-axis: Cumulative variance explained (dashed line + dots)
    geom_line(aes(y = CumulativeVariance * scale_factor, color = Species), linetype = "dashed", size = 1) +
    geom_point(aes(y = CumulativeVariance * scale_factor, color = Species), size = 2, shape = 21, fill = "white") +
    
    # Custom color mapping
    scale_color_manual(values = colors) +
    
    # Labels and theme
    labs(title = "Elbow Plot: Variance Explained by PCs",
         x = "Principal Component",
         y = "Percentage of Variance Explained") +
    
    # Secondary y-axis for cumulative VE (rescaled)
    scale_y_continuous(sec.axis = sec_axis(~ . / scale_factor, name = "Cumulative Variance Explained (%)")) +
    
    theme_minimal(base_size = 14)
  
  print(p)
}
