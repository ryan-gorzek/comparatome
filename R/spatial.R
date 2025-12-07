#' LoadSpatialData
#'
#' Load spatial transcriptomics data from 10x format or h5ad with coordinate information.
#' Creates a Seurat object with FOV (field of view) containing spatial centroids.
#'
#' @param data.path Character string specifying path to directory containing:
#'   - matrix.mtx.gz (or .mtx): Gene expression matrix
#'   - features.tsv.gz (or .tsv): Gene names/IDs  
#'   - barcodes.tsv.gz (or .tsv): Cell/spot barcodes
#'   - coords.tsv.gz (or .tsv): Spatial coordinates (barcodes, x, y columns)
#' @param project Character string for project name (default: "SpatialProject")
#' @param gene.column Integer specifying which column in features.tsv to use as gene names (default: 1)
#' @param assay Character string for assay name (default: "RNA")
#' @param fov.name Character string for FOV name (default: "FOV")
#'
#' @return Seurat object with:
#'   \itemize{
#'     \item RNA assay with count data
#'     \item Metadata columns: X, Y (spatial coordinates)
#'     \item FOV containing centroids for spatial visualization
#'   }
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Loads 10x format data using Seurat::Read10X
#'   \item Reads coordinate file (barcodes, x, y)
#'   \item Creates Seurat object with count matrix
#'   \item Adds X and Y coordinates to metadata
#'   \item Creates centroids and FOV for spatial plotting
#'   \item Subsets FOV to cells present in object
#' }
#'
#' Coordinate file format (coords.tsv):
#' ```
#' barcodes    x         y
#' BARCODE1    1234.5    5678.9
#' BARCODE2    1235.1    5679.2
#' ...
#' ```
#'
#' If coordinates need transformation (rotation, mirroring, etc.), apply these
#' before creating FOV or use TransformSpatialCoordinates().
#'
#' @export
#' @family spatial
#'
#' @examples
#' \dontrun{
#'   # Load spatial data from processed region
#'   obj <- LoadSpatialData(
#'     data.path = "path/to/region/COL1/",
#'     project = "V1_COL1",
#'     fov.name = "COL1"
#'   )
#'   
#'   # Verify FOV loaded
#'   ImageDimPlot(obj, fov = "COL1")
#' }
LoadSpatialData <- function(data.path, project = "SpatialProject", gene.column = 1, 
                           assay = "RNA", fov.name = "FOV") {
  
  library(Seurat)
  library(SeuratObject)
  
  # Load 10x format data
  obj.data <- Read10X(data.path, gene.column = gene.column)
  obj <- CreateSeuratObject(counts = obj.data, project = project, assay = assay)
  
  # Load coordinates
  coord.files <- list.files(data.path, pattern = "coords.*\\.tsv", full.names = TRUE)
  if (length(coord.files) == 0) {
    stop("Could not find coords.tsv file in ", data.path)
  }
  
  coords <- read.table(coord.files[1], header = TRUE, sep = "\t")
  
  # Match barcodes to object
  coords <- coords[match(colnames(obj), coords$barcodes), ]
  
  if (any(is.na(coords$x)) || any(is.na(coords$y))) {
    warning("Some cells missing coordinates - these will be removed")
    valid.cells <- !is.na(coords$x) & !is.na(coords$y)
    obj <- subset(obj, cells = colnames(obj)[valid.cells])
    coords <- coords[valid.cells, ]
  }
  
  # Add coordinates to metadata
  obj <- AddMetaData(obj, coords$x, "X")
  obj <- AddMetaData(obj, coords$y, "Y")
  
  # Create FOV
  cents.df <- data.frame(X = coords$y, Y = coords$x)
  rownames(cents.df) <- colnames(obj)
  cents <- CreateCentroids(cents.df)
  fov <- CreateFOV(
    cents,
    type = "centroids",
    assay = assay,
    key = Key(fov.name, quiet = TRUE)
  )
  
  # Subset to cells in object
  fov <- fov[Cells(obj)]
  obj[[fov.name]] <- fov
  
  return(obj)
}


#' TransformSpatialCoordinates
#'
#' Apply geometric transformations to spatial coordinates (rotation, scaling, translation).
#' Updates both metadata (X, Y columns) and FOV centroids.
#'
#' @param obj Seurat object with X and Y metadata columns
#' @param rotation_angle Numeric rotation angle in degrees (default: 0). 
#'   Positive values rotate counter-clockwise.
#' @param scale_x Numeric scaling factor for X axis (default: 1)
#' @param scale_y Numeric scaling factor for Y axis (default: 1)
#' @param translate_x Numeric translation offset for X axis (default: 0)
#' @param translate_y Numeric translation offset for Y axis (default: 0)
#' @param center_origin Logical whether to center coordinates at origin before transformation (default: FALSE)
#' @param fov.name Character string specifying which FOV to update (default: "FOV")
#'
#' @return Seurat object with transformed coordinates in:
#'   \itemize{
#'     \item Metadata: Updated X and Y columns
#'     \item FOV: Updated centroid coordinates
#'   }
#'   Original coordinates preserved as X_orig and Y_orig.
#'
#' @details
#' Transformation order:
#' \enumerate{
#'   \item Center at origin (if center_origin = TRUE)
#'   \item Apply rotation
#'   \item Apply scaling
#'   \item Apply translation
#' }
#'
#' Rotation matrix (counter-clockwise):
#' ```
#' [cos(θ)  -sin(θ)]
#' [sin(θ)   cos(θ)]
#' ```
#'
#' Common use cases:
#' - Align columns vertically: rotation_angle to make columnar axis vertical
#' - Normalize size: scale_x, scale_y to match reference dimensions
#' - Align multiple regions: translate_x, translate_y to position side-by-side
#'
#' @export
#' @family spatial
#'
#' @examples
#' \dontrun{
#'   # Rotate column by -7 degrees
#'   obj <- TransformSpatialCoordinates(
#'     obj, 
#'     rotation_angle = -7,
#'     fov.name = "COL1"
#'   )
#'   
#'   # Scale and align for multi-column comparison
#'   obj <- TransformSpatialCoordinates(
#'     obj,
#'     scale_x = 1.2,
#'     scale_y = 1.2, 
#'     translate_x = 1000,
#'     center_origin = TRUE
#'   )
#' }
TransformSpatialCoordinates <- function(obj, rotation_angle = 0, 
                                       scale_x = 1, scale_y = 1,
                                       translate_x = 0, translate_y = 0,
                                       center_origin = FALSE,
                                       fov.name = "FOV") {
  
  # Store original coordinates
  if (!"X_orig" %in% colnames(obj[[]])) {
    obj$X_orig <- obj$X
    obj$Y_orig <- obj$Y
  }
  
  # Get current coordinates
  coords <- as.matrix(obj[[]][, c("X", "Y")])
  
  # Center at origin if requested
  if (center_origin) {
    x_center <- mean(coords[, 1])
    y_center <- mean(coords[, 2])
    coords[, 1] <- coords[, 1] - x_center
    coords[, 2] <- coords[, 2] - y_center
  }
  
  # Apply rotation
  if (rotation_angle != 0) {
    theta <- rotation_angle * pi / 180
    rotation_matrix <- matrix(
      c(cos(theta), -sin(theta),
        sin(theta),  cos(theta)),
      nrow = 2, byrow = TRUE
    )
    coords <- coords %*% t(rotation_matrix)
  }
  
  # Apply scaling
  coords[, 1] <- coords[, 1] * scale_x
  coords[, 2] <- coords[, 2] * scale_y
  
  # Apply translation
  coords[, 1] <- coords[, 1] + translate_x
  coords[, 2] <- coords[, 2] + translate_y
  
  # Update metadata
  obj$X <- coords[, 1]
  obj$Y <- coords[, 2]
  
  # Update FOV if present
  if (fov.name %in% names(obj@images)) {
    cents.df <- data.frame(X = coords[, 2], Y = coords[, 1])
    rownames(cents.df) <- colnames(obj)
    cents <- CreateCentroids(cents.df)
    fov <- CreateFOV(
      cents,
      type = "centroids",
      assay = DefaultAssay(obj),
      key = Key(fov.name, quiet = TRUE)
    )
    fov <- fov[Cells(obj)]
    obj[[fov.name]] <- fov
  }
  
  return(obj)
}


#' AlignSpatialColumns
#'
#' Align multiple spatial regions (e.g., cortical columns) for side-by-side comparison.
#' Normalizes heights, scales to match reference, and positions horizontally with specified gaps.
#'
#' @param obj.list Named list of Seurat objects, each representing a spatial region
#' @param fov.names Character vector of FOV names corresponding to each object (default: names(obj.list))
#' @param reference_index Integer specifying which object to use as reference for height normalization (default: 1)
#' @param horizontal_gap Numeric gap between aligned regions in coordinate units (default: 0)
#' @param normalize_height Logical whether to scale all regions to match reference height (default: TRUE)
#' @param normalize_width Logical whether to scale all regions to match reference width (default: FALSE)
#'
#' @return List containing:
#'   \itemize{
#'     \item objects: List of transformed Seurat objects with updated coordinates
#'     \item combined: Merged Seurat object containing all regions
#'     \item metadata: Data frame with transformation parameters for each region
#'   }
#'
#' @details
#' Alignment workflow:
#' \enumerate{
#'   \item Calculate height (X range) and width (Y range) for each region
#'   \item If normalize_height, scale each region's Y coordinates to match reference height
#'   \item If normalize_width, scale each region's X coordinates to match reference width
#'   \item Shift all regions to start at Y = 0 (bottom alignment)
#'   \item Position regions horizontally with cumulative X offsets and gaps
#' }
#'
#' The combined object includes:
#' - All cells from all input objects
#' - Metadata column 'region' indicating source region
#' - Metadata column 'sample' if present in input objects
#' - Single merged FOV named "aligned"
#'
#' Useful for:
#' - Comparing cell type distributions across cortical columns
#' - Visualizing laminar organization side-by-side
#' - Density analysis across matched regions
#'
#' @export
#' @family spatial
#'
#' @examples
#' \dontrun{
#'   # Load multiple columns
#'   col1 <- LoadSpatialData("path/to/COL1/", fov.name = "COL1")
#'   col2 <- LoadSpatialData("path/to/COL2/", fov.name = "COL2")
#'   col3 <- LoadSpatialData("path/to/COL3/", fov.name = "COL3")
#'   
#'   # Align with normalized height and gaps
#'   aligned <- AlignSpatialColumns(
#'     obj.list = list(COL1 = col1, COL2 = col2, COL3 = col3),
#'     horizontal_gap = 500,
#'     normalize_height = TRUE
#'   )
#'   
#'   # Plot aligned columns
#'   ImageDimPlot(aligned$combined, fov = "aligned", group.by = "subclass")
#' }
AlignSpatialColumns <- function(obj.list, fov.names = NULL, reference_index = 1,
                               horizontal_gap = 0, normalize_height = TRUE, 
                               normalize_width = FALSE) {
  
  if (is.null(fov.names)) {
    fov.names <- names(obj.list)
  }
  
  if (length(obj.list) != length(fov.names)) {
    stop("Length of obj.list must equal length of fov.names")
  }
  
  # Calculate dimensions for each region
  dimensions <- lapply(obj.list, function(obj) {
    list(
      height = max(obj$X) - min(obj$X),
      width = max(obj$Y) - min(obj$Y),
      x_min = min(obj$X),
      y_min = min(obj$Y)
    )
  })
  
  # Reference dimensions
  ref_height <- dimensions[[reference_index]]$height
  ref_width <- dimensions[[reference_index]]$width
  
  # Transform each object
  transformed_objs <- list()
  current_x_offset <- 0
  transform_metadata <- data.frame()
  
  for (i in seq_along(obj.list)) {
    obj <- obj.list[[i]]
    dim <- dimensions[[i]]
    
    # Calculate scale factors
    scale_y <- if (normalize_height) ref_height / dim$height else 1
    scale_x <- if (normalize_width) ref_width / dim$width else 1
    
    # Calculate translation (bottom-align Y, cumulative X offset)
    translate_y <- -dim$y_min
    translate_x <- current_x_offset - dim$x_min * scale_x
    
    # Apply transformation
    obj <- TransformSpatialCoordinates(
      obj,
      scale_x = scale_x,
      scale_y = scale_y,
      translate_x = translate_x,
      translate_y = translate_y,
      fov.name = fov.names[i]
    )
    
    # Add region identifier
    obj$region <- names(obj.list)[i]
    
    # Store metadata
    meta <- data.frame(
      region = names(obj.list)[i],
      scale_x = scale_x,
      scale_y = scale_y,
      translate_x = translate_x,
      translate_y = translate_y,
      final_width = (max(obj$X) - min(obj$X)),
      final_height = (max(obj$Y) - min(obj$Y))
    )
    transform_metadata <- rbind(transform_metadata, meta)
    
    # Update offset for next region
    current_x_offset <- max(obj$X) + horizontal_gap
    
    transformed_objs[[names(obj.list)[i]]] <- obj
  }
  
  # Combine objects
  combined <- merge(transformed_objs[[1]], y = transformed_objs[-1])
  
  # Create combined FOV
  cents.df <- data.frame(X = combined$X, Y = combined$Y)
  rownames(cents.df) <- colnames(combined)
  cents <- CreateCentroids(cents.df)
  fov <- CreateFOV(
    cents,
    type = "centroids",
    assay = DefaultAssay(combined),
    key = Key("aligned", quiet = TRUE)
  )
  fov <- fov[Cells(combined)]
  combined[["aligned"]] <- fov
  
  return(list(
    objects = transformed_objs,
    combined = combined,
    metadata = transform_metadata
  ))
}


#' IntegrateAndLabelSpatial
#'
#' Integrate spatial transcriptomics data with reference dataset and assign labels via nearest neighbors.
#' Processes spatial data in chunks to handle large datasets efficiently, with automatic handling of small final chunks.
#'
#' @param obj.reference Seurat object with reference annotations (snRNA-seq or scRNA-seq)
#' @param obj.spatial Seurat object with spatial transcriptomics data
#' @param label.col Character string specifying metadata column in reference containing labels to transfer (e.g., "subclass", "type")
#' @param chunk_size Integer number of spatial cells to process per integration iteration (default: 1000)
#' @param integration_resolution Numeric resolution for clustering after integration (default: 0.5)
#' @param nn_fraction Numeric minimum fraction of neighbors required for label assignment (default: 0.25)
#' @param nn_neighbors Integer number of nearest neighbors to consider (default: 100)
#' @param subsample_reference Logical whether to subsample reference to match spatial chunk size (default: FALSE)
#' @param min_chunk_size Integer minimum size for final chunk; smaller chunks are merged with previous chunk (default: 200)
#'
#' @return Spatial Seurat object with new metadata columns:
#'   \itemize{
#'     \item [label.col]_nn: Nearest neighbor assigned labels
#'     \item [label.col]_batch: Batch number for each cell (tracking which integration chunk it was processed in)
#'   }
#'
#' @details
#' Integration workflow for large spatial datasets:
#' \enumerate{
#'   \item If total cells < chunk_size, process all cells in single batch
#'   \item Otherwise, split spatial object into chunks of size chunk_size
#'   \item If final chunk < min_chunk_size, merge with previous chunk to avoid processing too few cells
#'   \item For each chunk:
#'     a. Integrate with reference using SCTransform anchors
#'     b. Cluster integrated space at specified resolution  
#'     c. Assign labels via LabelByNearestNeighbors
#'     d. Extract assigned labels and batch numbers
#'   \item Combine labels from all chunks
#'   \item Add to original spatial object metadata
#' }
#'
#' Chunked processing advantages:
#' - Handles spatial datasets too large for single integration
#' - Reduces memory usage
#' - Provides progress tracking
#' - Automatically handles edge cases (small datasets, small final chunks)
#' - Tracks batch effects via batch column
#'
#' Label transfer uses nearest neighbors in integrated PCA space rather than
#' reference mapping, which works better for spatial data where technical
#' differences may be larger.
#'
#' The batch column allows checking for batch-specific biases in label assignment
#' and can be used as a covariate in downstream analyses.
#'
#' @export
#' @family spatial
#'
#' @examples
#' \dontrun{
#'   # Load reference with annotations
#'   obj.ref <- readRDS("opossum_v1_glutamatergic_processed.rds")
#'   
#'   # Load spatial data
#'   obj.spatial <- LoadSpatialData("path/to/spatial/region/")
#'   
#'   # Integrate and assign labels (standard)
#'   obj.spatial <- IntegrateAndLabelSpatial(
#'     obj.reference = obj.ref,
#'     obj.spatial = obj.spatial,
#'     label.col = "subclass",
#'     chunk_size = 1000,
#'     nn_fraction = 0.25,
#'     nn_neighbors = 100
#'   )
#'   
#'   # Check label distribution
#'   table(obj.spatial$subclass_nn)
#'   
#'   # Check batch distribution
#'   table(obj.spatial$subclass_batch)
#'   
#'   # Check for batch effects
#'   table(obj.spatial$subclass_nn, obj.spatial$subclass_batch)
#'   
#'   # Small dataset (processed in single batch)
#'   obj.small <- IntegrateAndLabelSpatial(
#'     obj.reference = obj.ref,
#'     obj.spatial = obj.spatial.small,  # <1000 cells
#'     label.col = "class",
#'     chunk_size = 1000
#'   )
#' }
IntegrateAndLabelSpatial <- function(obj.reference, obj.spatial, label.col = "subclass",
                                     chunk_size = 1000, integration_resolution = 0.5,
                                     nn_fraction = 0.25, nn_neighbors = 100,
                                     subsample_reference = FALSE,
                                     min_chunk_size = 200) {
  
  # Add method identifier
  obj.reference$method <- "reference"
  obj.spatial$method <- "spatial"
  
  # Initialize storage
  label_col_nn <- paste0(label.col, "_nn")
  label_df <- data.frame(
    cell = character(),
    label = character(),
    batch = integer(),
    stringsAsFactors = FALSE
  )
  
  # Get list of spatial cell names
  spatial_cells <- colnames(obj.spatial)
  n_cells <- length(spatial_cells)
  
  # If total cells < chunk_size, process all at once
  if (n_cells <= chunk_size) {
    cat(sprintf("Processing all %d cells in single batch\n", n_cells))
    n_chunks <- 1
  } else {
    n_chunks <- ceiling(n_cells / chunk_size)
    
    # Merge small final chunk with previous chunk
    final_chunk_size <- n_cells %% chunk_size
    if (final_chunk_size > 0 && final_chunk_size < min_chunk_size) {
      n_chunks <- n_chunks - 1
      cat(sprintf("Merging small final chunk (%d cells) with previous chunk\n", final_chunk_size))
    }
    
    cat(sprintf("Processing %d cells in %d chunks\n", n_cells, n_chunks))
  }
  
  # Process in chunks
  for (i in 1:n_chunks) {
    
    # Define chunk
    if (i == n_chunks) {
      # Last chunk gets all remaining cells
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- n_cells
    } else {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_cells)
    }
    
    chunk_cells <- spatial_cells[start_idx:end_idx]
    
    cat(sprintf("Chunk %d/%d: cells %d-%d (%d cells)\n", 
                i, n_chunks, start_idx, end_idx, length(chunk_cells)))
    
    # Subset spatial object
    obj.spatial.chunk <- subset(obj.spatial, cells = chunk_cells)
    
    # Optionally subsample reference
    if (subsample_reference) {
      ref_cells <- sample(colnames(obj.reference), min(ncol(obj.reference), chunk_size))
      obj.ref.chunk <- subset(obj.reference, cells = ref_cells)
    } else {
      obj.ref.chunk <- obj.reference
    }
    
    # Integrate
    obj.integrated <- IntegrateObjects(
      obj.ref.chunk, 
      obj.spatial.chunk,
      resolutions = integration_resolution,
      subsample = FALSE
    )
    
    # Assign labels via nearest neighbors
    obj.integrated <- LabelByNearestNeighbors(
      obj.integrated,
      ident = label.col,
      output_col = label_col_nn,
      fraction = nn_fraction,
      n.neighbors = nn_neighbors
    )
    
    # Extract labels for spatial cells
    spatial_labels <- obj.integrated[[label_col_nn]][obj.integrated$method == "spatial", , drop = FALSE]
    chunk_df <- data.frame(
      cell = rownames(spatial_labels),
      label = spatial_labels[, 1],
      batch = i,
      stringsAsFactors = FALSE
    )
    
    label_df <- rbind(label_df, chunk_df)
    
    cat(sprintf("  Assigned: %d cells, None: %.1f%%\n", 
                nrow(chunk_df),
                100 * sum(chunk_df$label == "None") / nrow(chunk_df)))
  }
  
  # Add labels and batch to original spatial object
  obj.spatial[[label_col_nn]] <- NA
  obj.spatial[[label_col_nn]][label_df$cell] <- label_df$label
  
  batch_col <- paste0(label.col, "_batch")
  obj.spatial[[batch_col]] <- NA
  obj.spatial[[batch_col]][label_df$cell] <- label_df$batch
  
  cat(sprintf("\nFinal label distribution:\n"))
  print(table(obj.spatial[[label_col_nn]], useNA = "ifany"))
  
  cat(sprintf("\nCells per batch:\n"))
  print(table(obj.spatial[[batch_col]]))
  
  return(obj.spatial)
}


#' ApplySpatialSelection
#'
#' Apply pre-defined spatial region selections from JSON parameter file to Seurat object.
#' Recreates selections made in Python interactive notebook.
#'
#' @param obj Seurat object with X and Y spatial coordinates in metadata
#' @param params_file Character string path to JSON file containing selection parameters
#' @param selection_names Character vector of specific selection names to apply (default: NULL for all)
#' @param selection_type Character string specifying type to extract: "regions", "sections", or "all" (default: "all")
#'
#' @return Named list of Seurat objects, one per selected region/section
#'
#' @details
#' JSON parameter file structure:
#' ```json
#' {
#'   "regions": [
#'     {
#'       "name": "COL1",
#'       "type": "rotated_rectangle",
#'       "angle_deg": -7,
#'       "x1": -7000, "x2": -5000,
#'       "y1": -14600, "y2": -13050
#'     }
#'   ],
#'   "sections": [
#'     {
#'       "name": "section1",
#'       "type": "polygon",
#'       "vertices": [[x1, y1], [x2, y2], ...]
#'     }
#'   ]
#' }
#' ```
#'
#' Selection types:
#' - **rotated_rectangle**: Applies rotation, then checks bounds
#' - **polygon**: Checks point-in-polygon using spatial coordinates
#'
#' Coordinates must match the coordinate system in obj$X and obj$Y.
#'
#' @export
#' @family spatial
#'
#' @examples
#' \dontrun{
#'   # Load full spatial dataset
#'   obj <- LoadSpatialData("path/to/spatial_raw/")
#'   
#'   # Apply selections from Python notebook
#'   regions <- ApplySpatialSelection(
#'     obj,
#'     params_file = "selection_params.json",
#'     selection_type = "regions"
#'   )
#'   
#'   # Access individual regions
#'   col1 <- regions$COL1
#'   col2 <- regions$COL2
#' }
ApplySpatialSelection <- function(obj, params_file, selection_names = NULL, 
                                 selection_type = "all") {
  
  library(jsonlite)
  library(sp)
  
  # Load parameters
  params <- fromJSON(params_file)
  
  # Determine which selections to process
  selections <- list()
  if (selection_type %in% c("regions", "all") && !is.null(params$regions)) {
    selections <- c(selections, params$regions)
  }
  if (selection_type %in% c("sections", "all") && !is.null(params$sections)) {
    selections <- c(selections, params$sections)
  }
  
  if (length(selections) == 0) {
    stop("No selections found in parameter file")
  }
  
  # Filter by name if specified
  if (!is.null(selection_names)) {
    selection_names_in_file <- sapply(selections, function(s) s$name)
    selections <- selections[selection_names_in_file %in% selection_names]
  }
  
  # Apply each selection
  results <- list()
  for (selection in selections) {
    
    cat(sprintf("Applying selection: %s (%s)\n", selection$name, selection$type))
    
    if (selection$type == "rotated_rectangle") {
      mask <- .get_rotated_rectangle_mask(obj, selection)
    } else if (selection$type == "polygon") {
      mask <- .get_polygon_mask(obj, selection)
    } else {
      warning(sprintf("Unknown selection type: %s", selection$type))
      next
    }
    
    if (sum(mask) > 0) {
      obj_subset <- subset(obj, cells = colnames(obj)[mask])
      obj_subset[[selection$name]] <- selection$name
      results[[selection$name]] <- obj_subset
      cat(sprintf("  Selected %d cells\n", sum(mask)))
    } else {
      warning(sprintf("Selection '%s' matched 0 cells", selection$name))
    }
  }
  
  return(results)
}

# Helper function for rotated rectangle mask
.get_rotated_rectangle_mask <- function(obj, selection) {
  # Get coordinates
  coords <- as.matrix(obj[[]][, c("X", "Y")])
  
  # Rotation matrix
  theta <- selection$angle_deg * pi / 180
  rotation_matrix <- matrix(
    c(cos(theta), -sin(theta),
      sin(theta),  cos(theta)),
    nrow = 2, byrow = TRUE
  )
  
  # Rotate coordinates
  rotated_coords <- coords %*% t(rotation_matrix)
  
  # Check bounds
  mask <- (
    rotated_coords[, 1] >= selection$x1 & 
    rotated_coords[, 1] <= selection$x2 &
    rotated_coords[, 2] >= selection$y1 & 
    rotated_coords[, 2] <= selection$y2
  )
  
  return(mask)
}

# Helper function for polygon mask
.get_polygon_mask <- function(obj, selection) {
  # Get coordinates  
  coords <- as.matrix(obj[[]][, c("X", "Y")])
  
  # Create polygon
  poly_coords <- do.call(rbind, selection$vertices)
  poly <- Polygon(poly_coords)
  poly <- Polygons(list(poly), ID = "1")
  poly <- SpatialPolygons(list(poly))
  
  # Check point in polygon
  points <- SpatialPoints(coords)
  mask <- !is.na(over(points, poly))
  
  return(mask)
}
