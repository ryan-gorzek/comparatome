#' MapObjects
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the mapping family.
#' @param seurat_obj1 (auto) parameter
#' @param seurat_obj2 (auto) parameter
#' @param idents (auto) parameter
#' @param assay (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family mapping
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
MapObjects <- function(seurat_obj1, seurat_obj2, idents, assay = "SCT") {
  
  objs <- list(seurat_obj1, seurat_obj2)
  # Perform SCTransform v2 on each object
  if (assay == "integrated") {
    objs <- lapply(objs, function(x) {
      # SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
      x <- RunPCA(x, npcs = 30, assay = "integrated", verbose = FALSE) %>%
           RunUMAP(reduction = "pca", dims = 1:30, assay = "integrated", return.model = TRUE, verbose = FALSE)
      return(x)
    })
  } else if (assay == "SCT") {
    objs <- lapply(objs, function(x) {
      x <- SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
           RunPCA(npcs = 30, verbose = FALSE) %>%
           RunUMAP(reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)
      return(x)
    })
  }
  objs.idx <- list(qry = c(1, 2), ref = c(2, 1))
  objs.mapped <- list()
  
  refdata <- list()
  for (id in idents) { refdata[[id]] <- id }
  
  for (idx in c(1, 2)) {

    # Transfer data from reference to query
    reference <- objs[[objs.idx$ref[idx]]]
    query <- objs[[objs.idx$qry[idx]]]
    anchors.query <- FindTransferAnchors(reference = reference, query = query, reference.reduction = "pca", dims = 1:30, k.filter = NA)
    # if (nrow(anchors.query@anchors) < 50) { k.weight = 10 } # floor(nrow(anchors.query@anchors) * 0.25)
    # else { k.weight = 50 }
    objs.mapped[[idx]] <- MapQuery(anchorset = anchors.query, 
                                   reference = reference, query = query, 
                                   refdata = refdata, 
                                   reference.reduction = "pca", 
                                   reduction.model = "umap")
  
  }
  return(objs.mapped)
}


#' MapObject
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the mapping family.
#' @param seurat_obj1 (auto) parameter
#' @param seurat_obj2 (auto) parameter
#' @param idents (auto) parameter
#' @param assay (auto) parameter
#' @param do.norm (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family mapping
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
MapObject <- function(seurat_obj1, seurat_obj2, idents, assay = "SCT", do.norm = TRUE) {
  
  objs <- list(seurat_obj1, seurat_obj2)
  # Perform SCTransform v2 on each object
  if (assay == "integrated") {
    objs <- lapply(objs, function(x) {
      # SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
      x <- RunPCA(x, npcs = 30, assay = "integrated", verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:30, assay = "integrated", return.model = TRUE, verbose = FALSE)
      return(x)
    })
  } else if (assay == "SCT") {
    objs <- lapply(objs, function(x) {
      x <- SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
        RunPCA(npcs = 30, verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)
      return(x)
    })
  } else if (do.norm == FALSE) { objs <- objs }
  
  refdata <- list()
  for (id in idents) { refdata[[id]] <- id }
    
  # Transfer data from reference to query
  reference <- objs[[1]]
  query <- objs[[2]]
  anchors.query <- FindTransferAnchors(reference = reference, query = query, reference.reduction = "pca", dims = 1:30)
  # if (nrow(anchors.query@anchors) < 50) { k.weight = 10 } # floor(nrow(anchors.query@anchors) * 0.25)
  # else { k.weight = 50 }
  obj.mapped <- MapQuery(anchorset = anchors.query, 
                         reference = reference, query = query, 
                         refdata = refdata, 
                         reference.reduction = "pca", 
                         reduction.model = "umap")
    
  return(obj.mapped)
}


#' PlotMapping
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the mapping family.
#' @param objs (auto) parameter
#' @param idents (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family mapping
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotMapping <- function(objs, idents = c("subclass", "type"), ident.order = NULL, title.key = "species") {
  
  for (obj.idx in c(1, 2)) {
    obj <- objs[[obj.idx]]
    for (id in idents) {
      id.levels <- as.character(unlist(unique(objs[[setdiff(c(1, 2), obj.idx)]][[id]])))
      for (red in c("umap", "ref.umap")) {
        for (cl in c("%s", "predicted.%s", "predicted.%s.score")) {
          if (cl == "predicted.%s.score") {
            plot <- FeaturePlot(obj, sprintf(cl, id), reduction = red) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
          }
          else {
            plot <- DimPlot(obj, reduction = red, group.by = sprintf(cl, id), label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
          }
          print(plot + ggtitle(paste(obj[[title.key]][1,], sprintf(cl, id), "on", red)))
        }
      }
      plot <- PlotMappedLabelsHeatmap(obj, id, id.levels, normalize = "row", ident.order = ident.order)
      print(plot)
      plot <- PlotMappingQualityHeatmap(obj, id, id.levels, sprintf("predicted.%s.score", id), ident.order = ident.order)
      print(plot)
    }
  }
  
}
