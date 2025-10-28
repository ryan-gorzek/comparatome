#' SubsampleObject
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the subsampling family.
#' @param seurat_obj (auto) parameter
#' @param metadata_col (auto) parameter
#' @param cells_per_category (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family subsampling
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
SubsampleObject <- function(seurat_obj, metadata_col, cells_per_category) {
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Get unique values in the specified metadata column
  unique_values <- unique(metadata[[metadata_col]])
  
  # Initialize a list to store subsampled cell names
  subsampled_cells <- list()
  
  # Loop through each unique value in the metadata column
  for (value in unique_values) {
    # Get cells that belong to the current metadata category
    cells_in_category <- rownames(metadata[metadata[[metadata_col]] == value, ])
    
    # Determine the number of cells to sample
    n_cells_to_sample <- min(length(cells_in_category), cells_per_category)
    
    # Sample cells
    sampled_cells <- sample(cells_in_category, n_cells_to_sample)
    
    # Add sampled cells to the list
    subsampled_cells <- c(subsampled_cells, sampled_cells)
  }
  
  # Subset the Seurat object to include only the subsampled cells
  subsampled_seurat_obj <- subset(seurat_obj, cells = as.character(subsampled_cells))
  
  return(subsampled_seurat_obj)
}


#' SubsampleObjectMultipleIterations
#'
#' Perform stratified subsampling of a Seurat object across multiple iterations, attempting to sample different cells in each iteration.
#' 
#' This function repeatedly subsamples cells from a Seurat object, stratified by a metadata column.
#' It attempts to minimize cell reuse across iterations by tracking previously sampled cells. Within each
#' iteration, the function samples `cells_per_category` cells from each unique value in the specified
#' metadata column. When a category lacks sufficient unsampled cells, the function draws from the
#' entire category pool to meet the target. This is useful for bootstrapping analyses, computing
#' confidence intervals, or assessing sampling variability in downstream analyses.
#' 
#' @param seurat_obj Seurat object to subsample.
#' @param metadata_col Character string specifying the name of a metadata column to stratify by (e.g., "cluster", "celltype").
#' @param cells_per_category Integer specifying the number of cells to sample from each unique value in `metadata_col` per iteration.
#' @param iterations Integer specifying the number of independent subsampling iterations to perform.
#' 
#' @details
#' The function implements a stratified subsampling strategy:
#' 
#' 1. **Within-iteration uniqueness**: Within each iteration, no cell is sampled more than once.
#' 
#' 2. **Cross-iteration diversity**: The function tracks cells sampled in previous iterations and attempts
#'    to select different cells in subsequent iterations. However, if a category has fewer total cells
#'    than `cells_per_category * iterations`, cells will necessarily be reused across iterations.
#' 
#' 3. **Category handling**: For each unique value in `metadata_col`, the function samples exactly
#'    `cells_per_category` cells per iteration. If a category contains fewer than `cells_per_category`
#'    cells, all available cells are sampled first, then additional cells are drawn from the entire
#'    category pool to reach the target.
#' 
#' 4. **Order preservation**: The order of categories matches the order encountered in the metadata column.
#' 
#' The returned list contains cell barcodes rather than Seurat objects to minimize memory usage, allowing
#' users to subset the original object as needed for each iteration.
#' 
#' @return List of length `iterations`, where each element is a character vector of cell barcodes
#'   representing the subsampled cells for that iteration. Each vector contains 
#'   `length(unique(metadata_col)) * cells_per_category` cell barcodes (assuming sufficient cells exist
#'   in each category).
#'   
#' @export
#' @family subsampling
#' 
#' @examples
#' \dontrun{
#'  # Perform 10 iterations of subsampling, taking 100 cells per cluster each time
#'  cell_lists <- SubsampleObjectMultipleIterations(
#'    seurat_obj = pbmc,
#'    metadata_col = "seurat_clusters",
#'    cells_per_category = 100,
#'    iterations = 10
#'  )
#'  
#'  # Subset the object for the first iteration
#'  pbmc_subsample_1 <- subset(pbmc, cells = cell_lists[[1]])
#'  
#'  # Use in a bootstrap workflow to compute standard errors
#'  marker_results <- lapply(cell_lists, function(cells) {
#'    obj_subsample <- subset(pbmc, cells = cells)
#'    FindMarkers(obj_subsample, ident.1 = "cluster1", ident.2 = "cluster2")
#'  })
#' }
SubsampleObjectMultipleIterations <- function(seurat_obj, metadata_col, cells_per_category, iterations) {
  
  # Extract metadata for cell lookup
  metadata <- seurat_obj@meta.data
  
  # Get unique values in the specified metadata column (e.g., cluster IDs, cell types)
  unique_values <- unique(metadata[[metadata_col]])
  
  # Initialize a list to store sampled cell barcodes for each iteration
  iterations_list <- vector("list", iterations)
  
  # Initialize a tracker to monitor which cells have been sampled in previous iterations
  # Structure: list of lists, where sampled_cells_overall[[iter]][[category]] contains cell barcodes
  sampled_cells_overall <- vector("list", iterations)
  
  # Iterate through each subsampling round
  for (iter in 1:iterations) {
    
    # Initialize a vector to accumulate sampled cells for this iteration
    subsampled_cells <- c()
    
    # Initialize a category-wise tracker for this iteration
    sampled_cells_iter <- list()
    
    # Loop through each category (e.g., each cluster or cell type)
    for (value in unique_values) {
      
      # Identify all cells belonging to the current category
      cells_in_category <- rownames(metadata[metadata[[metadata_col]] == value, ])
      
      # Extract cells previously sampled in prior iterations for this category (currently unused in logic below)
      previously_sampled <- unlist(lapply(sampled_cells_overall, function(x) if (!is.null(x[[value]])) x[[value]] else c()))
      
      # Compute cells available for sampling (exclude cells already sampled in this iteration)
      available_cells <- setdiff(cells_in_category, subsampled_cells)
      
      # Determine how many cells to sample
      if (length(available_cells) > 0) {
        
        if (length(available_cells) >= cells_per_category) {
          # Sufficient cells available: sample the target number without replacement
          sampled_cells <- sample(available_cells, cells_per_category, replace = FALSE)
          
        } else {
          # Insufficient cells available: take all available cells first
          sampled_cells <- available_cells
          
          # Calculate how many more cells are needed to reach the target
          remaining_needed <- cells_per_category - length(available_cells)
          
          # Define a pool of additional cells (all cells in category minus those just sampled)
          additional_pool <- setdiff(cells_in_category, sampled_cells)
          
          # Sample additional cells if pool is non-empty
          if (length(additional_pool) > 0) {
            additional_cells <- sample(additional_pool, remaining_needed, replace = FALSE)
            sampled_cells <- c(sampled_cells, additional_cells)
          }
        }
        
      } else {
        # No available cells (all have been sampled within this iteration already)
        # Fall back to sampling from the entire category without replacement
        sampled_cells <- sample(cells_in_category, cells_per_category, replace = FALSE)
      }
      
      # Accumulate sampled cells for this iteration
      subsampled_cells <- c(subsampled_cells, sampled_cells)
      
      # Track which cells were sampled from this category in this iteration
      sampled_cells_iter[[value]] <- sampled_cells
    }
    
    # Store the cell barcodes sampled in this iteration
    iterations_list[[iter]] <- subsampled_cells
    
    # Update the overall tracker with this iteration's sampling record
    sampled_cells_overall[[iter]] <- sampled_cells_iter
  }
  
  return(iterations_list)
}


#' SubsampleClasses
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the subsampling family.
#' @param seurat_obj (auto) parameter
#' @param meta_column (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family subsampling
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
SubsampleClasses <- function(seurat_obj, meta_column) {
  # Get the metadata
  metadata <- seurat_obj@meta.data
  
  # Determine the class sizes
  class_sizes <- table(metadata[[meta_column]])
  
  # Find the minimum class size
  min_size <- min(class_sizes)
  
  # Subsample each class to the minimum size
  sampled_cells <- unlist(lapply(names(class_sizes), function(class_name) {
    class_cells <- rownames(metadata[metadata[[meta_column]] == class_name, ])
    sample(class_cells, min_size)
  }))
  
  # Subset the Seurat object to the sampled cells
  downsampled_obj <- subset(seurat_obj, cells = sampled_cells)
  
  return(downsampled_obj)
}


#' SplitObjectHalf
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the subsampling family.
#' @param seurat_object (auto) parameter
#' @param seed (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family subsampling
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
SplitObjectHalf <- function(seurat_object, seed = 123) {
  # Get the number of cells
  num_cells <- ncol(seurat_object)
  
  # Generate a random split
  set.seed(seed)  # Setting seed for reproducibility
  random_split <- sample(1:num_cells, num_cells, replace = FALSE)
  
  # Split into two groups
  half1 <- random_split[1:(num_cells %/% 2)]
  half2 <- random_split[((num_cells %/% 2) + 1):num_cells]
  
  # Subset the Seurat object into two halves
  seurat_half1 <- subset(seurat_object, cells = colnames(seurat_object)[half1])
  seurat_half2 <- subset(seurat_object, cells = colnames(seurat_object)[half2])
  
  # Return a list containing the two halves
  return(list(obj1 = seurat_half1, obj2 = seurat_half2))
}


#' MatchDistribution
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the subsampling family.
#' @param query_obj (auto) parameter
#' @param reference_obj (auto) parameter
#' @param n_samples (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family subsampling
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
MatchDistribution <- function(query_obj, reference_obj, n_samples = NULL) {
  # Extract metadata
  query_data <- query_obj@meta.data %>% dplyr::select(nFeature_RNA, nCount_RNA)
  reference_data <- reference_obj@meta.data %>% dplyr::select(nFeature_RNA, nCount_RNA)
  
  # Determine the number of samples to draw if not specified
  if (is.null(n_samples)) {
    n_samples <- nrow(reference_data)
  }
  
  # Estimate the density of the reference dataset using 2D KDE
  kde <- kde2d(
    x = reference_data$nFeature_RNA, 
    y = reference_data$nCount_RNA, 
    n = 100
  )
  
  # Convert KDE results into a probability distribution
  kde_prob <- kde$z / sum(kde$z)
  
  # Sample from the KDE to get target distributions
  sampled_points <- MASS::mvrnorm(
    n = n_samples, 
    mu = c(mean(reference_data$nFeature_RNA), mean(reference_data$nCount_RNA)), 
    Sigma = cov(reference_data)
  )
  
  # Find the closest unique points in the query dataset to the sampled points without replacement
  available_indices <- seq_len(nrow(query_data))
  sampled_indices <- vector("integer", n_samples)
  
  for (i in 1:n_samples) {
    if (length(available_indices) == 0) break
    # Calculate distances between the sampled point and available query data
    distances <- sqrt(
      (query_data$nFeature_RNA[available_indices] - sampled_points[i, 1])^2 + 
        (query_data$nCount_RNA[available_indices] - sampled_points[i, 2])^2
    )
    # Find the closest available point
    closest_index <- available_indices[which.min(distances)]
    sampled_indices[i] <- closest_index
    # Remove the selected index from the available pool
    available_indices <- setdiff(available_indices, closest_index)
  }
  
  # Subset the query object to include only the sampled cells
  query_sampled <- subset(query_obj, cells = rownames(query_data)[sampled_indices])
  
  return(query_sampled)
}


#' ShuffleExpression
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the subsampling family.
#' @param seurat_obj (auto) parameter
#' @param metadata_column (auto) parameter
#' @param assay (auto) parameter
#' @param ignore_group (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family subsampling
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
ShuffleExpression <- function(seurat_obj, metadata_column, assay = "RNA", ignore_group = FALSE) {
  # Check if the metadata_column exists in the Seurat object metadata
  if (!(metadata_column %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", metadata_column, "not found in Seurat object metadata"))
  }
  
  # Create a copy of the expression data from the counts slot
  expr_data <- as.matrix(GetAssayData(seurat_obj, assay = assay, slot = "counts"))
  
  DefaultAssay(seurat_obj) <- assay
  Idents(seurat_obj) <- metadata_column
  
  if (ignore_group) {
    # Shuffle expression values across all cells
    shuffled_expr <- t(apply(expr_data, 1, function(x) x[sample(length(x))]))
    colnames(shuffled_expr) <- colnames(expr_data)
    
    # Create a new Seurat object to ensure the original remains unchanged
    seurat_obj_copy <- seurat_obj
    seurat_obj_copy <- SetAssayData(seurat_obj_copy, assay = assay, slot = "counts", new.data = as(shuffled_expr, "dgCMatrix"))
    
    return(seurat_obj_copy)
  } else {
    # Get unique categories
    categories <- unique(seurat_obj@meta.data[[metadata_column]])
    
    # Initialize a list to store shuffled expression matrices
    shuffled_expr_list <- list()
    cell_names_list <- list()
    
    # Loop through each category and shuffle expression values within that category
    for (category in categories) {
      print(category)
      
      # Get the cells belonging to the current category
      category_cells <- WhichCells(seurat_obj, idents = category)
      
      # Get the expression matrix for the current category cells
      category_expr <- expr_data[, category_cells]
      
      # Shuffle expression values for each gene within the current category
      shuffled_expr <- t(apply(category_expr, 1, function(x) x[sample(length(x))]))
      
      # Add the shuffled expression matrix to the list
      shuffled_expr_list[[category]] <- shuffled_expr
      cell_names_list[[category]] <- colnames(category_expr)
    }
    
    # Concatenate the shuffled expression matrices
    shuffled_expr_concat <- do.call(cbind, shuffled_expr_list)
    
    # Combine the cell names to maintain order
    all_cell_names <- unlist(cell_names_list)
    
    # Reorder the columns to match the original cell names
    cell_order <- colnames(expr_data)
    shuffled_expr_final <- shuffled_expr_concat[, match(cell_order, all_cell_names)]
    colnames(shuffled_expr_final) <- cell_order
    
    # Create a new Seurat object to ensure the original remains unchanged
    seurat_obj_copy <- seurat_obj
    seurat_obj_copy <- SetAssayData(seurat_obj_copy, assay = assay, slot = "counts", new.data = as(shuffled_expr_final, "dgCMatrix"))
    
    return(seurat_obj_copy)
  }
}
