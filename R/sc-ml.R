#' SaveIdentConfusionMatrices
#'
#' Generate and save confusion matrices for within-species classification using 50-50 
#' train-test splits. For each subclass and clustering resolution, cells are randomly 
#' split in half, with one half used as reference for MapQuery and the other half as 
#' query. The resulting confusion matrix shows classification accuracy when predicting 
#' cluster identity from transcriptomic profiles.
#'
#' @param obj Seurat object containing preprocessed single-cell data with clustering results
#' @param subclass.labels Character vector of subclass names to analyze
#' @param ident.labels Character vector of clustering resolution identifiers (e.g., "SCT_snn_res.0.2")
#' @param savepath Base directory path where confusion matrix plots will be saved
#' @param assay Assay to use for mapping (default: "SCT")
#' @param n_iters Number of train-test iterations to average (default: 1)
#'
#' @return NULL (saves plots to disk)
#'
#' @details
#' The function iterates through each subclass and clustering resolution:
#' \itemize{
#'   \item Extracts cells belonging to the current subclass
#'   \item Randomly splits cells 50-50 into reference and query sets
#'   \item Uses MapObject to train on reference and predict on query
#'   \item Generates confusion matrix via PlotMappedLabelsHeatmap
#'   \item Saves plot to [savepath]/[subclass]/[resolution]/[subclass.col]_Mapping.png
#' }
#'
#' Output directory structure: savepath/SubclassName/ResolutionName/
#'
#' @export
#' @family ml-confusion
#'
#' @examples
#' \dontrun{
#'   SaveIdentConfusionMatrices(
#'     obj = seurat_obj,
#'     subclass.labels = c("IT_A", "IT_B", "L5PT"),
#'     ident.labels = c("SCT_snn_res.0.2", "SCT_snn_res.0.5"),
#'     savepath = "output/confusion_matrices/"
#'   )
#' }
SaveIdentConfusionMatrices <- function(obj, subclass.labels, ident.labels, savepath, assay = "SCT", n_iters = 1) {
  
  DefaultAssay(obj) <- assay
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, gsub("/", "", sbcl), "/")
    make_folder(folder_path)
    
    for (id in ident.labels) {
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj) <- sbcl.col
      if (sbcl %in% levels(obj)) {
        
        obj.sbcl.id <- subset(obj, idents = sbcl)
        
        DefaultAssay(obj.sbcl.id) <- assay
        Idents(obj.sbcl.id) <- id
        levels(obj.sbcl.id) <- sort_idents(levels(obj.sbcl.id))
        if (length(levels(obj.sbcl.id)) > 1) {
          
          id.path <- paste0(folder_path, gsub("/", "", id), "/")
          make_folder(id.path)
          
          # Perform 50-50 train-test split and mapping
          set.seed(42)
          all_cells <- Cells(obj.sbcl.id)
          train_cells <- sample(all_cells, size = floor(length(all_cells) / 2))
          test_cells <- setdiff(all_cells, train_cells)
          
          obj.train <- subset(obj.sbcl.id, cells = train_cells)
          obj.test <- subset(obj.sbcl.id, cells = test_cells)
          
          # Map test to train
          obj.mapped <- MapObject(obj.train, obj.test, idents = id, assay = assay, do.norm = FALSE)
          
          # Generate confusion matrix
          id.levels <- levels(Idents(obj.sbcl.id))
          confusion_plot <- PlotMappedLabelsHeatmap(obj.mapped, id, id.levels, normalize = "row")
          
          if (!is.null(confusion_plot)) {
            confusion_plot <- confusion_plot +
              coord_equal() + 
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_rect(fill = "white", color = NA), 
                    plot.background = element_rect(fill = "white", color = NA))
            
            ggsave(paste0(id.path, sbcl.col, "_Mapping.png"), plot = confusion_plot, width = 5, height = 5, dpi = 300)
          }
        }
      }
    }
  }
}


#' SaveSubclassConfusionMatrices
#'
#' Generate and save confusion matrices for within-species subclass classification using 
#' 50-50 train-test splits. Similar to SaveIdentConfusionMatrices but operates at the 
#' subclass level rather than on individual clustering resolutions.
#'
#' @param obj Seurat object containing preprocessed single-cell data with subclass annotations
#' @param subclass.cols Character, name of metadata column containing subclass labels
#' @param subclass.order Character vector specifying the order of subclasses for plotting
#' @param savepath Base directory path where confusion matrix plots will be saved
#' @param return.plot Logical, whether to return the plot object (default: FALSE)
#' @param colormap_upper_limit Numeric, upper limit for color scale (default: NULL)
#' @param assay Assay to use for mapping (default: "SCT")
#'
#' @return NULL or ggplot object if return.plot = TRUE
#'
#' @details
#' Creates confusion matrix showing classification accuracy across subclasses using
#' 50-50 train-test split with MapObject workflow. Useful for assessing overall
#' subclass separability and classification robustness.
#'
#' @export
#' @family ml-confusion
#'
#' @examples
#' \dontrun{
#'   SaveSubclassConfusionMatrices(
#'     obj = seurat_obj,
#'     subclass.cols = "subclass",
#'     subclass.order = c("IT_A", "IT_B", "L5PT", "L6CT"),
#'     savepath = "output/confusion_matrices/"
#'   )
#' }
SaveSubclassConfusionMatrices <- function(obj, subclass.cols, subclass.order, savepath, return.plot = FALSE, colormap_upper_limit = NULL, assay = "SCT") {
  
  DefaultAssay(obj) <- assay
  
  for (sbcl in subclass.cols) {
    
    folder_path <- paste0(savepath)
    make_folder(folder_path)
    
    DefaultAssay(obj) <- assay
    Idents(obj) <- sbcl
    levels(obj) <- sort_by_reference(levels(obj), subclass.order)
    
    # Perform 50-50 train-test split and mapping
    set.seed(42)
    all_cells <- Cells(obj)
    train_cells <- sample(all_cells, size = floor(length(all_cells) / 2))
    test_cells <- setdiff(all_cells, train_cells)
    
    obj.train <- subset(obj, cells = train_cells)
    obj.test <- subset(obj, cells = test_cells)
    
    # Map test to train
    obj.mapped <- MapObject(obj.train, obj.test, idents = sbcl, assay = assay, do.norm = FALSE)
    
    # Generate confusion matrix
    sbcl.levels <- levels(Idents(obj))
    confusion_plot <- PlotMappedLabelsHeatmap(obj.mapped, sbcl, sbcl.levels, normalize = "row")
    
    if (!is.null(confusion_plot)) {
      confusion_plot <- confusion_plot +
        coord_equal() + 
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_rect(fill = "white", color = NA), 
              plot.background = element_rect(fill = "white", color = NA))
      
      ggsave(paste0(folder_path, sbcl, "_Mapping.png"), plot = confusion_plot, width = 5, height = 5, dpi = 300)
      
      if (return.plot) { return(confusion_plot) }
    }
  }
}
