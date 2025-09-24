#' SaveIdentConfusionMatrices
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the ml-confusion family.
#' @param obj (auto) parameter
#' @param subclass.labels (auto) parameter
#' @param ident.labels (auto) parameter
#' @param savepath (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family ml-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
SaveIdentConfusionMatrices <- function(obj, subclass.labels, ident.labels, savepath) {
  
  DefaultAssay(obj) <- "SCT"
  
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
        
        DefaultAssay(obj.sbcl.id) <- "SCT"
        Idents(obj.sbcl.id) <- id
        levels(obj.sbcl.id) <- sort_idents(levels(obj.sbcl.id))
        if (length(levels(obj.sbcl.id)) > 1) {
          
          id.path <- paste0(folder_path, gsub("/", "", id), "/")
          make_folder(id.path)

          mdl <- TrainModel(obj.sbcl.id, training_genes = VariableFeatures(obj))
          
          if (!is.null(mdl$confusion)) {
            confusion_plot <- mdl$confusion +
              coord_equal() + 
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_rect(fill = "white", color = NA), 
                    plot.background = element_rect(fill = "white", color = NA))
            
            ggsave(paste0(id.path, sbcl.col, "_XGBoost.png"), plot = confusion_plot, width = 5, height = 5, dpi = 300)
          }
        }
      }
    }
  }
}


#' SaveSubclassConfusionMatrices
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the ml-confusion family.
#' @param obj (auto) parameter
#' @param subclass.cols (auto) parameter
#' @param subclass.order (auto) parameter
#' @param savepath (auto) parameter
#' @param return.plot (auto) parameter
#' @param colormap_upper_limit (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family ml-confusion
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
SaveSubclassConfusionMatrices <- function(obj, subclass.cols, subclass.order, savepath, return.plot = FALSE, colormap_upper_limit = NULL) {
  
  DefaultAssay(obj) <- "SCT"
  
  for (sbcl in subclass.cols) {
    
    folder_path <- paste0(savepath)
    make_folder(folder_path)
    
    DefaultAssay(obj) <- "SCT"
    Idents(obj) <- sbcl
    levels(obj) <- sort_by_reference(levels(obj), subclass.order)
    
    mdl <- TrainModel(obj, training_genes = VariableFeatures(obj))
    
    if (!is.null(mdl$confusion)) {
      confusion_plot <- mdl$confusion +
        coord_equal() + 
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_rect(fill = "white", color = NA), 
              plot.background = element_rect(fill = "white", color = NA))
      
      if (!is.null(colormap_upper_limit)) {
        confusion_plot <- confusion_plot +
          scale_fill_gradient(limits = c(0, colormap_upper_limit), 
                              low = "white", 
                              high = "red", 
                              oob = scales::squish)
      } else {
        confusion_plot <- confusion_plot +
          scale_fill_gradient(low = "white", high = "red")
      }
      
      ggsave(paste0(folder_path, sbcl, "_XGBoost.png"), plot = confusion_plot, width = 5, height = 5, dpi = 300)
      if (return.plot == TRUE) {
        return(confusion_plot)
      }
    }
  }
}
