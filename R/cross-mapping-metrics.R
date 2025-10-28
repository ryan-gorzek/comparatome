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
#' @param col.low Color for low values in heatmap gradient (default: "white")
#' @param col.high Color for high values in heatmap gradient (default: "red")
#' @param x.lab.rot Logical, whether to rotate x-axis labels 90 degrees (default: TRUE)
#'
#' @return ggplot2 object showing confusion matrix heatmap with percentage/count labels in tiles
#'
#' @details
#' The function expects the Seurat object to contain:
#' \itemize{
#'   \item Original labels in metadata column specified by column_name
#'   \item Predicted labels in column named "predicted.[column_name]" (created by MapObject)
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
#'   print(p)
#' }
PlotMappedLabelsHeatmap <- function(data, column_name, column_levels, normalize = NULL, ident.order = NULL, col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  
  # Create confusion matrix
  confusion_matrix <- table(as.character(unlist(data[[column_name]])), as.character(unlist(data[[paste0("predicted.", column_name)]])))
  confusion_matrix <- as.matrix(confusion_matrix)
  
  add_zeros_to_table <- function(tbl, new_row_names, new_col_names) {
    # Add new rows of zeros
    for (row_name in new_row_names) {
      if (!any(row_name %in% rownames(tbl))) {
        tbl <- rbind(tbl, setNames(t(rep(0, ncol(tbl))), row_name))
        orig_names <- rownames(tbl)[rownames(tbl) != ""]
        rownames(tbl) <- c(orig_names, row_name)
      }
    }
    
    # Add new columns of zeros
    for (col_name in new_col_names) {
      if (!any(col_name %in% colnames(tbl))) {
        tbl <- cbind(tbl, setNames(rep(0, nrow(tbl)), col_name))
        orig_names <- colnames(tbl)[colnames(tbl) != ""]
        colnames(tbl) <- c(orig_names, col_name)
      }
    }
    
    return(tbl)
  }
  
  # Ensure all possible levels are present in the confusion matrix
  row_levels <- as.character(unlist(unique(data[[column_name]])))
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
  
  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(unlist(df$ident))
  }
  
  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)
  
  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_ident(row_levels, unique(melted$row))
    col_levels <- sort_ident(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_ident(row_levels, ident.order)
    col_levels <- sort_ident(col_levels, ident.order)
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


#' PlotSubsampledMappedLabelsHeatmap
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the mapping-metrics family.
#' @param true.labels (auto) parameter
#' @param pred.labels (auto) parameter
#' @param column_levels (auto) parameter
#' @param normalize (auto) parameter
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
PlotSubsampledMappedLabelsHeatmap <- function(true.labels, pred.labels, column_levels, normalize = NULL, ident.order = NULL, col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  
  # Create confusion matrix
  confusion_matrix <- table(as.character(unlist(true.labels)), as.character(unlist(pred.labels)))
  confusion_matrix <- as.matrix(confusion_matrix)
  
  add_zeros_to_table <- function(tbl, new_row_names, new_col_names) {
    # Add new rows of zeros
    for (row_name in new_row_names) {
      if (!any(row_name %in% rownames(tbl))) {
        tbl <- rbind(tbl, setNames(t(rep(0, ncol(tbl))), row_name))
        orig_names <- rownames(tbl)[rownames(tbl) != ""]
        rownames(tbl) <- c(orig_names, row_name)
      }
    }
    
    # Add new columns of zeros
    for (col_name in new_col_names) {
      if (!any(col_name %in% colnames(tbl))) {
        tbl <- cbind(tbl, setNames(rep(0, nrow(tbl)), col_name))
        orig_names <- colnames(tbl)[colnames(tbl) != ""]
        colnames(tbl) <- c(orig_names, col_name)
      }
    }
    
    return(tbl)
  }
  
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
    if (normalize == "row") {
      melted <- ddply(melted, .(row), transform, Percentage = Count / sum(Count) * 100)
    } else if (normalize == "col") {
      melted <- ddply(melted, .(col), transform, Percentage = Count / sum(Count) * 100)
    }
  } else {
    melted$Percentage <- melted$Count
  }
  
  # # Sorting function
  # sort_ident <- function(ident, primary_order) {
  #   primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
  #   suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
  #   suffix_numeric <- suppressWarnings(as.numeric(suffix))
  #   suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
  #   suffix_numeric[is.na(suffix_numeric)] <- Inf
  #   df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
  #   df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
  #   return(unlist(df$ident))
  # }
  
  sort_ident <- function(ident, primary_order) {
    # Ensure that all identifiers in 'ident' are part of 'primary_order'
    ident <- factor(ident, levels = primary_order, ordered = TRUE)
    return(as.character(levels(ident)))  # Return the factor levels in the specified order
  }
  
  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)

  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_ident(row_levels, unique(melted$row))
    col_levels <- sort_ident(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_ident(row_levels, ident.order)
    col_levels <- sort_ident(col_levels, ident.order)
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
  average_matrix <<- with(data, tapply(as.numeric(unlist(data[[value_column]])), 
                                       list(as.character(unlist(data[[column_name]])), 
                                       as.character(unlist(data[[paste0("predicted.", column_name)]]))), mean, na.rm = TRUE))
  average_matrix[is.na(average_matrix)] <- 0
  average_matrix <- average_matrix * 100

  add_zeros_to_table <- function(tbl, new_row_names, new_col_names) {
    # Add new rows of zeros
    for (row_name in new_row_names) {
      if (!any(row_name %in% rownames(tbl))) {
        tbl <- rbind(tbl, setNames(t(rep(0, ncol(tbl))), row_name))
        orig_names <- rownames(tbl)[rownames(tbl) != ""]
        rownames(tbl) <- c(orig_names, row_name)
      }
    }

    # Add new columns of zeros
    for (col_name in new_col_names) {
      if (!any(col_name %in% colnames(tbl))) {
        tbl <- cbind(tbl, setNames(rep(0, nrow(tbl)), col_name))
        orig_names <- colnames(tbl)[colnames(tbl) != ""]
        colnames(tbl) <- c(orig_names, col_name)
      }
    }

    return(tbl)
  }

  # Ensure all possible levels are present in the average matrix
  row_levels <<- as.character(unlist(unique(data[[column_name]])))
  col_levels <<- column_levels
  rows_to_add <<- row_levels[row_levels %in% rownames(average_matrix) == FALSE]
  cols_to_add <<- col_levels[col_levels %in% colnames(average_matrix) == FALSE]
  average_matrix <- add_zeros_to_table(average_matrix, rows_to_add, cols_to_add)

  # average_df <- as.data.frame(average_matrix)

  # Melt the dataframe for ggplot2
  melted <- melt(average_matrix)
  colnames(melted) <- c("row", "col", "Value")
  melted$Value <- as.numeric(melted$Value)

  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(unlist(df$ident))
  }

  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)

  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_ident(row_levels, unique(melted$row))
    col_levels <- sort_ident(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_ident(row_levels, ident.order)
    col_levels <- sort_ident(col_levels, ident.order)
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
  
  # Function to add zeros to the confusion matrix if levels are missing
  add_zeros_to_table <- function(tbl, new_row_names, new_col_names) {
    for (row_name in new_row_names) {
      if (!any(row_name %in% rownames(tbl))) {
        tbl <- rbind(tbl, setNames(t(rep(0, ncol(tbl))), row_name))
        orig_names <- rownames(tbl)[rownames(tbl) != ""]
        rownames(tbl) <- c(orig_names, row_name)
      }
    }
    
    for (col_name in new_col_names) {
      if (!any(col_name %in% colnames(tbl))) {
        tbl <- cbind(tbl, setNames(rep(0, nrow(tbl)), col_name))
        orig_names <- colnames(tbl)[colnames(tbl) != ""]
        colnames(tbl) <- c(orig_names, col_name)
      }
    }
    
    return(tbl)
  }
  
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
