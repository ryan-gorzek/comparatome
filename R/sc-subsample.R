#' SubsampleObject
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the subsampling family.
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
#' Auto-generated roxygen skeleton for comparatome.
Part of the subsampling family.
#' @param seurat_obj (auto) parameter
#' @param metadata_col (auto) parameter
#' @param cells_per_category (auto) parameter
#' @param iterations (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family subsampling
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
SubsampleObjectMultipleIterations <- function(seurat_obj, metadata_col, cells_per_category, iterations) {
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Get unique values in the specified metadata column
  unique_values <- unique(metadata[[metadata_col]])
  
  # Initialize a list to store sampled cells for each iteration
  iterations_list <- vector("list", iterations)
  
  # Initialize a list to keep track of already sampled cells
  sampled_cells_overall <- vector("list", iterations)
  
  for (iter in 1:iterations) {
    # Initialize a vector to store sampled cells for this iteration
    subsampled_cells <- c()
    
    # Initialize a list to track sampled cells for this iteration
    sampled_cells_iter <- list()
    
    for (value in unique_values) {
      # Get cells that belong to the current metadata category
      cells_in_category <- rownames(metadata[metadata[[metadata_col]] == value, ])
      
      # Exclude cells that have already been sampled in previous iterations
      previously_sampled <- unlist(lapply(sampled_cells_overall, function(x) if (!is.null(x[[value]])) x[[value]] else c()))
      available_cells <- setdiff(cells_in_category, subsampled_cells) # Ensure no duplicates within iteration
      
      # Check if there are available cells to sample
      if (length(available_cells) > 0) {
        if (length(available_cells) >= cells_per_category) {
          # Sample without replacement
          sampled_cells <- sample(available_cells, cells_per_category, replace = FALSE)
        } else {
          # Sample all available cells
          sampled_cells <- available_cells
          
          # If needed, sample additional cells from the entire pool
          remaining_needed <- cells_per_category - length(available_cells)
          
          # Sample from the previously used cells plus the remaining available cells
          additional_pool <- setdiff(cells_in_category, sampled_cells)
          
          # Check if additional_pool is non-empty before sampling
          if (length(additional_pool) > 0) {
            additional_cells <- sample(additional_pool, remaining_needed, replace = FALSE)
            # Combine the available and additional sampled cells
            sampled_cells <- c(sampled_cells, additional_cells)
          }
        }
      } else {
        # If no available cells, sample with replacement from the entire category
        sampled_cells <- sample(cells_in_category, cells_per_category, replace = FALSE)
      }
      
      # Add sampled cells to the subsampled_cells vector
      subsampled_cells <- c(subsampled_cells, sampled_cells)
      
      # Store sampled cells for this category in the tracker for this iteration
      sampled_cells_iter[[value]] <- sampled_cells
    }
    
    # Store the subsampled cells for this iteration
    iterations_list[[iter]] <- subsampled_cells
    
    # Store the sampled cells for this iteration in the overall tracker
    sampled_cells_overall[[iter]] <- sampled_cells_iter
  }
  
  return(iterations_list)
}


#' SubsampleClasses
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the subsampling family.
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
Part of the subsampling family.
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
Part of the subsampling family.
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
Part of the subsampling family.
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
