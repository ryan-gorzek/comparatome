#' PlotMappedLabelsHeatmap
#'
#' Generate confusion matrix heatmap comparing true labels to predicted labels from cell mapping.
#' Creates a tile plot showing classification accuracy with optional row or column normalization.
#' Used internally by SaveIdentConfusionMatrices and SaveSubclassConfusionMatrices.
#'
#' @param data Seurat object containing both original labels and predicted labels from MapObject
#' @param column_name Character, name of metadata column containing true labels
#' @param column_levels Character vector specifying all possible label levels (ensures complete matrix even if some labels absent)
#' @param normalize Character, normalization method: "row" (percentage within true label), "col" (percentage within predicted label), 
#'   or NULL for raw counts (default: NULL)
#' @param ident.order Character vector specifying custom order for labels in plot. If NULL, uses alphanumeric sorting with 
#'   intelligent handling of suffixes (default: NULL)
#' @param true.column Character, optional name of metadata column containing true labels. If NULL, uses column_name (default: NULL)
#' @param pred.column Character, optional name of metadata column containing predicted labels. If NULL, uses "predicted.[column_name]" (default: NULL)
#' @param col.low Color for low values in heatmap gradient (default: "white")
#' @param col.high Color for high values in heatmap gradient (default: "red")
#' @param x.lab.rot Logical, whether to rotate x-axis labels 90 degrees (default: TRUE)
#'
#' @return ggplot2 object showing confusion matrix heatmap with percentage/count labels in tiles
#'
#' @details
#' The function expects the Seurat object to contain:
#' \itemize{
#'   \item Original labels in metadata column specified by column_name or true.column
#'   \item Predicted labels in column named "predicted.[column_name]" (created by MapObject) or pred.column
#' }
#'
#' Matrix construction:
#' \enumerate{
#'   \item Creates confusion matrix from true vs predicted labels
#'   \item Adds zero-filled rows/columns for any missing levels from column_levels
#'   \item Normalizes by row or column if specified
#'   \item Sorts labels intelligently (primary label, then numeric suffix, then alphabetic)
#' }
#'
#' Text size automatically adjusts based on matrix dimensions (smaller for >10 classes).
#' Plot uses fixed aspect ratio for square tiles.
#'
#' @keywords internal
#' @family mapping-metrics
#'
#' @examples
#' \dontrun{
#'   # After mapping with MapObject
#'   obj.mapped <- MapObject(obj.train, obj.test, idents = "subclass")
#'   
#'   # Generate confusion matrix
#'   p <- PlotMappedLabelsHeatmap(
#'     data = obj.mapped,
#'     column_name = "subclass",
#'     column_levels = c("IT_A", "IT_B", "L5PT", "L6CT"),
#'     normalize = "row"
#'   )
#'   
#'   # With custom column names
#'   p <- PlotMappedLabelsHeatmap(
#'     data = obj.mapped,
#'     column_name = "cluster",
#'     column_levels = c("IT_A", "IT_B", "L5PT", "L6CT"),
#'     true.column = "SCT_snn_res.0.2",
#'     pred.column = "predicted.subclass",
#'     normalize = "row"
#'   )
#'   print(p)
#' }
PlotMappedLabelsHeatmap <- function(data, column_name, column_levels, normalize = NULL, ident.order = NULL, true.column = NULL, pred.column = NULL, col.low = "white", col.high = "red", x.lab.rot = TRUE) {

  # Use custom column names if provided, otherwise default to column_name and predicted.column_name
  true.col <- if (!is.null(true.column)) true.column else column_name
  pred.col <- if (!is.null(pred.column)) pred.column else paste0("predicted.", column_name)
  
  # Create confusion matrix
  confusion_matrix <- table(as.character(unlist(data[[true.col]])), as.character(unlist(data[[pred.col]])))
  confusion_matrix <- as.matrix(confusion_matrix)
  
  # Ensure all possible levels are present in the confusion matrix
  row_levels <- as.character(unlist(unique(data[[true.col]])))
  col_levels <- column_levels
  rows_to_add <- row_levels[row_levels %in% rownames(confusion_matrix) == FALSE]
  cols_to_add <- col_levels[col_levels %in% colnames(confusion_matrix) == FALSE]
  confusion_matrix <- add_zeros_to_table(confusion_matrix, rows_to_add, cols_to_add)
  
  confusion_df <- as.data.frame(confusion_matrix)
  
  # Melt the dataframe for ggplot2
  melted <- melt(confusion_matrix)
  colnames(melted) <- c("row", "col", "Count")
  melted$Count <- as.numeric(melted$Count)
  
  # Normalize if needed
  if (!is.null(normalize)) {
    if (normalize %in% c("row", "col")) {
      melted <- ddply(melted, normalize, transform, Percentage = Count / sum(Count) * 100)
    } else {
      melted$Percentage <- melted$Count
    }
  }
  
  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)
  
  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_idents(row_levels, unique(melted$row))
    col_levels <- sort_idents(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_idents(row_levels, ident.order)
    col_levels <- sort_idents(col_levels, ident.order)
  }
  
  melted$row <- factor(melted$row, levels = row_levels)
  melted$col <- factor(melted$col, levels = col_levels)
  
  if (nrow(confusion_matrix) < 10) { fontsize = 5 }
  else { fontsize = 4 }
  
  # Plot heatmap
  p <- ggplot(melted, aes(y = factor(row, levels = rev(row_levels)), x = factor(col, levels = col_levels))) + 
    geom_tile(aes(fill = Percentage)) + 
    scale_fill_gradient(low = col.low, high = col.high, limits = c(0, max(melted$Percentage))) + 
    geom_text(aes(label = sprintf("%.1f", Percentage)), size = fontsize) +
    theme_bw() + 
    ylab(true.col) + 
    xlab(pred.col) + 
    theme(
      axis.text.x = element_text(size = 16, face = "italic", hjust = 0, angle = ifelse(x.lab.rot, 90, 0)),
      axis.text.y = element_text(size = 16, face = "italic"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    ) +
    coord_fixed()
  
  # Ensure correct rotation
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  }
  
  return(p)
}


#' PlotSubsampledMappedLabelsHeatmap
#'
#' Generate confusion matrix heatmap from aggregated subsampled mapping results.
#' Similar to PlotMappedLabelsHeatmap but operates on pre-aggregated label vectors 
#' rather than a Seurat object, allowing visualization of results from multiple mapping iterations.
#'
#' @param true.labels Character vector of true labels (e.g., cluster IDs from query dataset)
#' @param pred.labels Character vector of predicted labels (e.g., subclass assignments from reference dataset)
#'   Must be same length as true.labels
#' @param column_levels Character vector specifying all possible predicted label levels to include in matrix,
#'   ensures complete matrix even if some labels absent from predictions
#' @param normalize Character, normalization method: "row" (percentage within true label), 
#'   "col" (percentage within predicted label), or NULL for raw counts (default: NULL)
#' @param ident.order Numeric/character vector specifying custom display order for both rows and columns.
#'   If NULL, uses factor-based ordering from column_levels (default: NULL)
#' @param col.low Color for low values in heatmap gradient (default: "white")
#' @param col.high Color for high values in heatmap gradient (default: "red")
#' @param x.lab.rot Logical, whether to rotate x-axis labels 90 degrees (default: TRUE)
#'
#' @return ggplot2 object showing confusion matrix heatmap with percentage/count labels in tiles
#'
#' @details
#' Designed for visualizing cross-species mapping results from multiple subsampled iterations.
#' Typical workflow:
#' \enumerate{
#'   \item Subsample query and reference datasets multiple times
#'   \item Map query to reference for each iteration using MapObject
#'   \item Aggregate all true and predicted labels across iterations
#'   \item Visualize with this function to show overall mapping patterns
#' }
#'
#' Matrix construction:
#' \enumerate{
#'   \item Creates confusion matrix from true vs predicted label vectors
#'   \item Adds zero-filled rows/columns for missing levels from column_levels
#'   \item Normalizes by row or column if specified
#'   \item Orders labels according to ident.order (or uses factor levels if not specified)
#' }
#'
#' Text size automatically adjusts (smaller for >10 classes). Fixed aspect ratio ensures square tiles.
#' Unlike PlotMappedLabelsHeatmap, this function doesn't extract labels from a Seurat object,
#' making it suitable for pre-aggregated results stored as vectors.
#'
#' @export
#' @family mapping-metrics
#'
#' @examples
#' \dontrun{
#'   # After multiple subsampled mapping iterations
#'   mapping.results <- list(
#'     true.labels = c(),
#'     pred.labels = c()
#'   )
#'   
#'   for (i in 1:10) {
#'     obj.query.sub <- SubsampleObject(obj.query, "cluster", 100)
#'     obj.ref.sub <- SubsampleObject(obj.ref, "subclass", 100)
#'     obj.mapped <- MapObject(obj.ref.sub, obj.query.sub, idents = "subclass")
#'     
#'     mapping.results$true.labels <- c(
#'       mapping.results$true.labels,
#'       as.character(obj.mapped$cluster)
#'     )
#'     mapping.results$pred.labels <- c(
#'       mapping.results$pred.labels,
#'       as.character(obj.mapped$predicted.subclass)
#'     )
#'   }
#'   
#'   # Visualize aggregated results
#'   p <- PlotSubsampledMappedLabelsHeatmap(
#'     true.labels = mapping.results$true.labels,
#'     pred.labels = mapping.results$pred.labels,
#'     column_levels = c("IT_A", "IT_B", "L5PT", "L6CT"),
#'     normalize = "row",
#'     ident.order = c(1:20, "IT_A", "IT_B", "L5PT", "L6CT")
#'   )
#'   print(p)
#' }
PlotSubsampledMappedLabelsHeatmap <- function(true.labels, pred.labels, column_levels, normalize = NULL, ident.order = NULL, col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  
  # Create confusion matrix
  confusion_matrix <- table(as.character(unlist(true.labels)), as.character(unlist(pred.labels)))
  confusion_matrix <- as.matrix(confusion_matrix)
  
  # Ensure all possible levels are present in the confusion matrix
  row_levels <- as.character(unlist(unique(true.labels)))
  col_levels <- column_levels
  rows_to_add <- row_levels[row_levels %in% rownames(confusion_matrix) == FALSE]
  cols_to_add <- col_levels[col_levels %in% colnames(confusion_matrix) == FALSE]
  confusion_matrix <- add_zeros_to_table(confusion_matrix, rows_to_add, cols_to_add)
  
  confusion_df <- as.data.frame(confusion_matrix)
  
  # Melt the dataframe for ggplot2
  melted <- melt(confusion_matrix)
  colnames(melted) <- c("row", "col", "Count")
  melted$Count <- as.numeric(melted$Count)
  
  # Normalize if needed
  if (!is.null(normalize)) {
    if (normalize %in% c("row", "col")) {
      melted <- ddply(melted, normalize, transform, Percentage = Count / sum(Count) * 100)
    } else {
      melted$Percentage <- melted$Count
    }
  }
  
  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)

  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_idents(row_levels, unique(melted$row))
    col_levels <- sort_idents(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_idents(row_levels, ident.order)
    col_levels <- sort_idents(col_levels, ident.order)
  }
  
  melted$row <- factor(melted$row, levels = row_levels)
  melted$col <- factor(melted$col, levels = col_levels)
  
  if (nrow(confusion_matrix) < 10) { fontsize = 5 }
  else { fontsize = 4 }
  
  # Plot heatmap
  p <- ggplot(melted, aes(y = factor(row, levels = rev(row_levels)), x = factor(col, levels = col_levels))) + 
    geom_tile(aes(fill = Percentage)) + 
    scale_fill_gradient(low = col.low, high = col.high, limits = c(0, max(melted$Percentage))) + 
    geom_text(aes(label = sprintf("%.1f", Percentage)), size = fontsize) +
    theme_bw() + 
    theme(
      axis.text.x = element_text(size = 16, face = "italic", hjust = 0, angle = ifelse(x.lab.rot, 90, 0)),
      axis.text.y = element_text(size = 16, face = "italic"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    ) +
    coord_fixed()
  
  # Ensure correct rotation
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  }
  
  return(p)
}


#' PlotMappingQualityHeatmap
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the mapping-metrics family.
#' @param data (auto) parameter
#' @param column_name (auto) parameter
#' @param column_levels (auto) parameter
#' @param value_column (auto) parameter
#' @param ident.order (auto) parameter
#' @param col.low (auto) parameter
#' @param col.high (auto) parameter
#' @param x.lab.rot (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family mapping-metrics
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotMappingQualityHeatmap <- function(data, column_name, column_levels, value_column, ident.order = NULL, col.low = "white", col.high = "red", x.lab.rot = TRUE) {

  # Create average matrix
  average_matrix <- with(data, tapply(as.numeric(unlist(data[[value_column]])), 
                                       list(as.character(unlist(data[[column_name]])), 
                                       as.character(unlist(data[[paste0("predicted.", column_name)]]))), mean, na.rm = TRUE))
  average_matrix[is.na(average_matrix)] <- 0
  average_matrix <- average_matrix * 100

  # Ensure all possible levels are present in the average matrix
  row_levels <- as.character(unlist(unique(data[[column_name]])))
  col_levels <- column_levels
  rows_to_add <- row_levels[row_levels %in% rownames(average_matrix) == FALSE]
  cols_to_add <- col_levels[col_levels %in% colnames(average_matrix) == FALSE]
  average_matrix <- add_zeros_to_table(average_matrix, rows_to_add, cols_to_add)

  # average_df <- as.data.frame(average_matrix)

  # Melt the dataframe for ggplot2
  melted <- melt(average_matrix)
  colnames(melted) <- c("row", "col", "Value")
  melted$Value <- as.numeric(melted$Value)

  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)

  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_idents(row_levels, unique(melted$row))
    col_levels <- sort_idents(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_idents(row_levels, ident.order)
    col_levels <- sort_idents(col_levels, ident.order)
  }

  melted$row <- factor(melted$row, levels = row_levels)
  melted$col <- factor(melted$col, levels = col_levels)

  if (nrow(average_matrix) < 10) { fontsize = 5 }
  else { fontsize = 4 }

  # Plot heatmap
  p <- ggplot(melted, aes(y = factor(row, levels = rev(row_levels)), x = factor(col, levels = col_levels))) +
    geom_tile(aes(fill = Value)) +
    scale_fill_gradient(low = col.low, high = col.high, limits = c(0, max(melted$Value, na.rm = TRUE))) +
    geom_text(aes(label = ifelse(Value == 0, "NA", sprintf("%.0f", Value))), size = fontsize) +
    theme_bw() +
    ylab(column_name) +
    xlab(paste0("predicted_", column_name)) +
    theme(
      axis.text.x = element_text(size = 16, face = "italic", hjust = 0, angle = ifelse(x.lab.rot, 90, 0)),
      axis.text.y = element_text(size = 16, face = "italic"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    ) +
    coord_fixed()

  # Ensure correct rotation
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  }

  return(p)
}


#' MappingAccuracy
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the mapping-metrics family.
#' @param data (auto) parameter
#' @param column_name (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family mapping-metrics
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
MappingAccuracy <- function(data, column_name) {
  
  # Create confusion matrix
  confusion_matrix <- table(as.character(unlist(data[[column_name]])), as.character(unlist(data[[paste0("predicted.", column_name)]])))
  confusion_matrix <- as.matrix(confusion_matrix)
  
  # Ensure all possible levels are present in the confusion matrix
  row_levels <- as.character(unlist(unique(data[[column_name]])))
  rows_to_add <- row_levels[row_levels %in% rownames(confusion_matrix) == FALSE]
  cols_to_add <- row_levels[row_levels %in% colnames(confusion_matrix) == FALSE]  # Using row_levels since column_levels is removed
  confusion_matrix <- add_zeros_to_table(confusion_matrix, rows_to_add, cols_to_add)
  
  # Calculate accuracy for each subclass
  subclass_accuracy <- sapply(rownames(confusion_matrix), function(subclass) {
    true_positives <- confusion_matrix[subclass, subclass]
    total_actual <- sum(confusion_matrix[subclass, ])
    if (total_actual > 0) {
      return(true_positives / total_actual)
    } else {
      return(NA)  # In case there are no actual instances of the subclass
    }
  })
  
  # Convert to dataframe
  accuracy_df <- data.frame(Subclass = names(subclass_accuracy), Accuracy = subclass_accuracy)
  
  # Return dataframe
  return(accuracy_df)
}
