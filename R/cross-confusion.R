#' PlotSubclassCrossConfusionMatrices
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the cross-confusion family.
#' @param obj1 (auto) parameter
#' @param obj2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param assay (auto) parameter
#' @param subclass.order (auto) parameter
#' @param n_iters (auto) parameter
#' @param all.genes (auto) parameter
#' @param upsample (auto) parameter
#' @param downsample (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family cross-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotSubclassCrossConfusionMatrices <- function(obj1, obj2, sample.name1, sample.name2, assay = "SCT", subclass.order, n_iters = 10, all.genes = FALSE, upsample = FALSE, downsample = FALSE) {
  library(gridExtra)
  
  DefaultAssay(obj1) <- assay
  DefaultAssay(obj2) <- assay
  Idents(obj1) <- "subclass"
  Idents(obj2) <- "subclass"
  levels(obj1) <- sort_by_reference(levels(obj1), subclass.order)
  levels(obj2) <- sort_by_reference(levels(obj2), subclass.order)
  objs.list <- list(obj1, obj2)
  
  if (all.genes) {
    variable.features = intersect(rownames(obj1), rownames(obj2))
  }
  else {
    variable.features <<- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
  }
  
  # Helper function to compute confusion matrix
  compute_confusion_matrix <- function(train_obj, test_obj) {
    mdl <- TrainModel(train_obj, assay = assay, training_genes = variable.features, upsample = upsample, downsample = downsample)
    confusion_matrix <- BuildConfusionMatrix(test_obj, train_obj, model = mdl, assay = assay)
    as.table(confusion_matrix)
  }
  
  # Initialize lists to store all row and column names
  all_row_names <- unique(c(levels(Idents(obj1)), levels(Idents(obj2))))
  all_col_names <- unique(c(levels(Idents(obj1)), levels(Idents(obj2))))
  
  # Function to align matrices to have the same row and column names
  align_matrix <- function(mat, row_names, col_names) {
    aligned_mat <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                          dimnames = list(row_names, col_names))
    mat_rownames <- rownames(mat)
    mat_colnames <- colnames(mat)
    for (r in mat_rownames) {
      for (c in mat_colnames) {
        aligned_mat[r, c] <- mat[r, c]
      }
    }
    return(aligned_mat)
  }
  
  # Function to filter rows and columns with non-zero sums
  filter_non_zero <- function(mat) {
    row_sums <- rowSums(mat)
    col_sums <- colSums(mat)
    non_zero_rows <- rownames(mat)[row_sums != 0]
    non_zero_cols <- colnames(mat)[col_sums != 0]
    mat[non_zero_rows, non_zero_cols]
  }
  
  # Initialize matrices to store summed percentages and lists to store individual matrices
  confusion_matrix1_sum <- matrix(0, nrow = length(all_row_names), ncol = length(all_col_names),
                                  dimnames = list(all_row_names, all_col_names))
  confusion_matrix2_sum <- matrix(0, nrow = length(all_row_names), ncol = length(all_col_names),
                                  dimnames = list(all_row_names, all_col_names))
  
  confusion_matrices1 <- list()
  confusion_matrices2 <- list()
  
  # Compute confusion matrices for both directions and store them
  for (i in 1:n_iters) {
    cm1 <- compute_confusion_matrix(objs.list[[1]], objs.list[[2]])
    cm2 <- compute_confusion_matrix(objs.list[[2]], objs.list[[1]])
    cm1_aligned <- align_matrix(cm1, all_row_names, all_col_names)
    cm2_aligned <- align_matrix(cm2, all_row_names, all_col_names)
    confusion_matrix1_sum <- confusion_matrix1_sum + cm1_aligned
    confusion_matrix2_sum <- confusion_matrix2_sum + cm2_aligned
    confusion_matrices1[[i]] <- cm1_aligned
    confusion_matrices2[[i]] <- cm2_aligned
  }
  
  # Average the confusion matrices
  confusion_matrix1_avg <- confusion_matrix1_sum / n_iters
  confusion_matrix2_avg <- confusion_matrix2_sum / n_iters
  
  # Filter rows and columns with non-zero sums
  confusion_matrix1_avg <- filter_non_zero(confusion_matrix1_avg)
  confusion_matrix2_avg <- filter_non_zero(confusion_matrix2_avg)
  
  # Create the averaged confusion matrix plots
  plot1 <- create_confusion_matrix_plot(confusion_matrix1_avg, sample.name1, sample.name2, subclass.order)
  plot1 <- plot1 + ggtitle(paste0(length(variable.features), " Features"))
  plot2 <- create_confusion_matrix_plot(confusion_matrix2_avg, sample.name2, sample.name1, subclass.order)
  plot2 <- plot2 + ggtitle(paste0(length(variable.features), " Features"))
  
  # Create the individual confusion matrix plots
  grid_plot1 <- plot_individual_confusion_matrices(confusion_matrices1, n_iters, sample.name1, sample.name2, subclass.order, filter_non_zero)
  grid_plot2 <- plot_individual_confusion_matrices(confusion_matrices2, n_iters, sample.name2, sample.name1, subclass.order, filter_non_zero)
  
  # Return the plots in a list
  return(list(avg_confusion_plot1 = plot1, avg_confusion_plot2 = plot2, individual_grid_plot1 = grid_plot1, individual_grid_plot2 = grid_plot2))
}


#' PlotIdentCrossConfusionMatrices
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the cross-confusion family.
#' @param obj1 (auto) parameter
#' @param obj2 (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param assay (auto) parameter
#' @param subclass.labels (auto) parameter
#' @param ident.labels (auto) parameter
#' @param n_iters (auto) parameter
#' @param all.genes (auto) parameter
#' @param ident.genes (auto) parameter
#' @param upsample (auto) parameter
#' @param downsample (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family cross-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotIdentCrossConfusionMatrices <- function(obj1, obj2, sample.name1, sample.name2, assay = "SCT", subclass.labels, ident.labels, n_iters = 10, all.genes = FALSE, ident.genes = FALSE, upsample = FALSE, downsample = FALSE) {
  library(gridExtra)
  
  DefaultAssay(obj1) <- assay
  DefaultAssay(obj2) <- assay
  
  objs.list <- list(obj1, obj2)
  if (all.genes) {
    variable.features = intersect(rownames(obj1), rownames(obj2))
  }
  else {
    variable.features <<- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
  }
  
  all_plots <- list()
  
  for (sbcl in subclass.labels) {
    
    all_plots[[sbcl]] <- list()
    
    for (id in ident.labels) {
      
      all_plots[[sbcl]][[id]] <- list(avg1 = NA, avg2 = NA, grid1 = NA, grid2 = NA)
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj1) <- sbcl.col
      Idents(obj2) <- sbcl.col
      
      if (sbcl %in% levels(obj1) & sbcl %in% levels(obj2)) {
        
        obj1.sbcl.id <- subset(obj1, idents = sbcl)
        obj2.sbcl.id <- subset(obj2, idents = sbcl)
        
        if (ident.genes == TRUE && all.genes == FALSE && assay == "SCT") {
          print("Using IntegrationFeatures from the Ident level...")
          obj1.sbcl.id.temp <- SCTransform(obj1.sbcl.id, verbose = FALSE)
          obj2.sbcl.id.temp <- SCTransform(obj2.sbcl.id, verbose = FALSE)
          objs.list <- list(obj1.sbcl.id.temp, obj2.sbcl.id.temp)
          variable.features <- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
        }
        else if (ident.genes == TRUE && all.genes == FALSE && assay == "integrated") {
          print("Using IntegrationFeatures from the Ident level...")
          obj1.sbcl.id.temp <- FindVariableFeatures(obj1.sbcl.id, assay = "integrated", nfeatures = 3000)
          obj2.sbcl.id.temp <- FindVariableFeatures(obj2.sbcl.id, assay = "integrated", nfeatures = 3000)
          objs.list <- list(obj1.sbcl.id.temp, obj2.sbcl.id.temp)
          variable.features <- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
        }
        
        DefaultAssay(obj1.sbcl.id) <- assay
        DefaultAssay(obj2.sbcl.id) <- assay
        
        Idents(obj1.sbcl.id) <- id
        Idents(obj2.sbcl.id) <- id
        
        levels(obj1.sbcl.id) <- sort_idents(levels(obj1.sbcl.id))
        levels(obj2.sbcl.id) <- sort_idents(levels(obj2.sbcl.id))
        
        if (length(levels(obj1.sbcl.id)) > 1 & length(levels(obj2.sbcl.id)) > 1) {
          
          confusion_matrices1 <- list()
          confusion_matrices2 <- list()

          names1 <- unique(levels(Idents(obj1.sbcl.id)))
          names2 <- unique(levels(Idents(obj2.sbcl.id)))
          
          # Initialize matrices to store summed percentages
          confusion_matrix1_sum <- matrix(0, nrow = length(names2), ncol = length(names1),
                                          dimnames = list(names2, names1))
          confusion_matrix2_sum <- matrix(0, nrow = length(names1), ncol = length(names2),
                                          dimnames = list(names1, names2))
          
          valid_matrices1 <- 0
          valid_matrices2 <- 0
          
          for (i in 1:n_iters) {
            
            # Helper function to compute confusion matrix
            compute_confusion_matrix <- function(train_obj, test_obj) {
              mdl <- TrainModel(train_obj, assay = assay, training_genes = variable.features, upsample = upsample, downsample = downsample)
              confusion_matrix <- BuildConfusionMatrix(test_obj, train_obj, model = mdl, assay = assay)
              as.table(confusion_matrix)
            }
            
            confusion_matrix1 <- compute_confusion_matrix(obj1.sbcl.id, obj2.sbcl.id)
            confusion_matrix2 <- compute_confusion_matrix(obj2.sbcl.id, obj1.sbcl.id)
            
            if (!is.null(confusion_matrix1)) {
              cm1 <- as.matrix(confusion_matrix1)
              if (nrow(cm1) > 0 && ncol(cm1) > 0) {
                cm1_aligned <- align_matrix(cm1, names2, names1)
                confusion_matrix1_sum <- confusion_matrix1_sum + cm1_aligned
                confusion_matrices1[[valid_matrices1 + 1]] <- cm1_aligned
                valid_matrices1 <- valid_matrices1 + 1
              }
            }
            
            if (!is.null(confusion_matrix2)) {
              cm2 <- as.matrix(confusion_matrix2)
              if (nrow(cm2) > 0 && ncol(cm2) > 0) {
                cm2_aligned <- align_matrix(cm2, names1, names2)
                confusion_matrix2_sum <- confusion_matrix2_sum + cm2_aligned
                confusion_matrices2[[valid_matrices2 + 1]] <- cm2_aligned
                valid_matrices2 <- valid_matrices2 + 1
              }
            }
          }
          
          # Average the confusion matrices if there are valid matrices
          if (valid_matrices1 > 0) {
            confusion_matrix1_avg <- confusion_matrix1_sum / valid_matrices1
            avg_plot1 <- create_ident_confusion_matrix_plot(confusion_matrix1_avg, sample.name1, sample.name2)
            avg_plot1 <- avg_plot1 + ggtitle(paste0(length(variable.features), " Features"))
            grid_plot1 <- plot_ident_individual_confusion_matrices(confusion_matrices1, valid_matrices1, sample.name1, sample.name2)
            all_plots[[sbcl]][[id]][["avg1"]] <- avg_plot1
            all_plots[[sbcl]][[id]][["grid1"]] <- grid_plot1
          }
          
          if (valid_matrices2 > 0) {
            confusion_matrix2_avg <- confusion_matrix2_sum / valid_matrices2
            avg_plot2 <- create_ident_confusion_matrix_plot(confusion_matrix2_avg, sample.name2, sample.name1)
            avg_plot2 <- avg_plot2 + ggtitle(paste0(length(variable.features), " Features"))
            grid_plot2 <- plot_ident_individual_confusion_matrices(confusion_matrices2, valid_matrices2, sample.name2, sample.name1)
            all_plots[[sbcl]][[id]][["avg2"]] <- avg_plot2
            all_plots[[sbcl]][[id]][["grid2"]] <- grid_plot2
          }
        }
      }
    }
  }
  
  return(all_plots)
}


#' create_ident_confusion_matrix_plot
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the cross-confusion family.
#' @param confusion_matrix (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family cross-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
create_ident_confusion_matrix_plot <- function(confusion_matrix, sample.name1, sample.name2) {
  row.levels <- sort_idents(rownames(confusion_matrix))
  col.levels <- sort_idents(colnames(confusion_matrix))
  
  if (!is.null(confusion_matrix)) {
    p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                             row.levels = row.levels, col.levels = col.levels)
    confusion_plot <- p +
      coord_equal() +
      labs(
        x = paste0("Predicted (Train = ", sample.name1, ")"),
        y = paste0("True (Test = ", sample.name2, ")")
      ) +
      theme(axis.text.x = element_text(angle = 90),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))
    return(confusion_plot)
  }
}


#' plot_ident_individual_confusion_matrices
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the cross-confusion family.
#' @param confusion_matrices (auto) parameter
#' @param n_valid_iters (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family cross-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
plot_ident_individual_confusion_matrices <- function(confusion_matrices, n_valid_iters, sample.name1, sample.name2) {
  plot_list <- list()
  
  for (i in 1:n_valid_iters) {
    confusion_matrix <- confusion_matrices[[i]]
    row.levels <- sort_idents(rownames(confusion_matrix))
    col.levels <- sort_idents(colnames(confusion_matrix))
    
    if (!is.null(confusion_matrix)) {
      p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                               row.levels = row.levels, col.levels = col.levels)
      confusion_plot <- p +
        coord_equal() +
        labs(
          title = paste0("Iteration ", i)
        ) +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      plot_list[[i]] <- confusion_plot
    }
  }
  
  grid_plot <- grid.arrange(grobs = plot_list, ncol = 5)
  return(grid_plot)
}


#' create_confusion_matrix_plot
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the cross-confusion family.
#' @param confusion_matrix (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family cross-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
create_confusion_matrix_plot <- function(confusion_matrix, sample.name1, sample.name2, subclass.order) {
  row.levels <- sort_by_reference(rownames(confusion_matrix), subclass.order)
  col.levels <- sort_by_reference(colnames(confusion_matrix), subclass.order)
  
  if (!is.null(confusion_matrix)) {
    p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                             row.levels = row.levels, col.levels = col.levels)
    confusion_plot <- p +
      coord_equal() +
      labs(
        x = paste0("Predicted (Train = ", sample.name1, ")"),
        y = paste0("True (Test = ", sample.name2, ")")
      ) +
      theme(axis.text.x = element_text(angle = 90),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))
    return(confusion_plot)
  }
}


#' plot_individual_confusion_matrices
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the cross-confusion family.
#' @param confusion_matrices (auto) parameter
#' @param n_iters (auto) parameter
#' @param sample.name1 (auto) parameter
#' @param sample.name2 (auto) parameter
#' @param subclass.order (auto) parameter
#' @param filter_func (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family cross-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
plot_individual_confusion_matrices <- function(confusion_matrices, n_iters, sample.name1, sample.name2, subclass.order, filter_func) {
  plot_list <- list()
  
  for (i in 1:n_iters) {
    confusion_matrix <- filter_func(confusion_matrices[[i]])
    row.levels <- sort_by_reference(rownames(confusion_matrix), subclass.order)
    col.levels <- sort_by_reference(colnames(confusion_matrix), subclass.order)
    
    if (!is.null(confusion_matrix)) {
      p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                               row.levels = row.levels, col.levels = col.levels)
      confusion_plot <- p +
        coord_equal() +
        labs(
          title = paste0("Iteration ", i)
        ) +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      plot_list[[i]] <- confusion_plot
    }
  }
  
  grid_plot <- grid.arrange(grobs = plot_list, ncol = 5)
  return(grid_plot)
}
