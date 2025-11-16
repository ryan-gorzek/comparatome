#' MapObjects
#'
#' Perform bidirectional label transfer between two Seurat objects using anchor-based mapping.
#' Maps labels in both directions simultaneously and returns both mapped query objects for
#' reciprocal comparison.
#'
#' @param seurat_obj1 First Seurat object with cell type labels in metadata
#' @param seurat_obj2 Second Seurat object with cell type labels in metadata
#' @param idents Character vector of metadata column names containing labels to transfer 
#'   bidirectionally (e.g., c("subclass", "type"))
#' @param assay Character, assay for anchor finding and mapping: "SCT" runs SCTransform v2 
#'   (default), "integrated" uses existing integrated assay and runs PCA/UMAP
#'
#' @return List of length 2 containing mapped query objects:
#' \itemize{
#'   \item [[1]]: seurat_obj2 with labels from seurat_obj1 (obj1→obj2 mapping)
#'   \item [[2]]: seurat_obj1 with labels from seurat_obj2 (obj2→obj1 mapping)
#' }
#' Each returned object includes:
#' \itemize{
#'   \item predicted.[ident]: Transferred labels for each ident column
#'   \item predicted.[ident].score: Confidence scores (0-1) for predictions
#'   \item ref.umap: Query cells projected into reference UMAP space
#' }
#'
#' @details
#' **Workflow for each mapping direction:**
#' \enumerate{
#'   \item Normalize both objects (SCTransform v2 or PCA on integrated assay)
#'   \item Run PCA (30 components) and UMAP with return.model = TRUE
#'   \item Find transfer anchors using reference PCA space (k.filter = NA for stability)
#'   \item Transfer labels from reference to query via MapQuery
#'   \item Project query cells into reference UMAP coordinates
#' }
#' 
#' **Normalization behavior:**
#' Always performs normalization. For pre-normalized objects, use MapObject with do.norm = FALSE.
#' 
#' **Use cases:**
#' - Cross-species reciprocal mapping (e.g., mouse↔opossum)
#' - Comparing mapping quality and symmetry in both directions
#' - Generating paired confusion matrices for bidirectional validation
#' - Assessing which species provides better reference for the other
#' 
#' **Key differences from MapObject:**
#' - **MapObjects**: Bidirectional (returns 2 objects), always normalizes, uses k.filter = NA
#' - **MapObject**: Unidirectional (returns 1 object), optional normalization, default k.filter
#' 
#' **Technical notes:**
#' - k.filter = NA prevents anchor filtering, useful when cell counts differ between objects
#' - UMAP return.model = TRUE enables projection of query into reference space
#' - Both objects should have similar feature sets for optimal anchor finding
#' - Memory intensive - considers running MapObject twice if memory limited
#'
#' @export
#' @family mapping
#'
#' @examples
#' \dontrun{
#'   # Bidirectional cross-species subclass mapping
#'   mapped <- MapObjects(
#'     seurat_obj1 = obj.mouse,
#'     seurat_obj2 = obj.opossum,
#'     idents = c("subclass"),
#'     assay = "SCT"
#'   )
#'   
#'   # Mouse labels transferred to opossum
#'   obj.opossum.mapped <- mapped[[1]]
#'   table(obj.opossum.mapped$predicted.subclass)
#'   summary(obj.opossum.mapped$predicted.subclass.score)
#'   
#'   # Opossum labels transferred to mouse
#'   obj.mouse.mapped <- mapped[[2]]
#'   table(obj.mouse.mapped$predicted.subclass)
#'   
#'   # Compare mapping symmetry
#'   DimPlot(mapped[[1]], reduction = "ref.umap", group.by = "predicted.subclass")
#'   DimPlot(mapped[[2]], reduction = "ref.umap", group.by = "predicted.subclass")
#'   
#'   # Visualize both mappings with helper function
#'   PlotMapping(mapped, idents = "subclass")
#' }
MapObjects <- function(seurat_obj1, seurat_obj2, idents, assay = "SCT") {
  
  library(Seurat)
  
  # Normalize and prepare both objects
  objs <- list(seurat_obj1, seurat_obj2)
  
  if (assay == "integrated") {
    objs <- lapply(objs, function(obj) {
      obj <- RunPCA(obj, npcs = 30, assay = "integrated", verbose = FALSE) %>%
             RunUMAP(reduction = "pca", dims = 1:30, assay = "integrated", 
                     return.model = TRUE, verbose = FALSE)
      return(obj)
    })
  } else if (assay == "SCT") {
    objs <- lapply(objs, function(obj) {
      obj <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = FALSE, 
                        verbose = FALSE) %>%
             RunPCA(npcs = 30, verbose = FALSE) %>%
             RunUMAP(reduction = "pca", dims = 1:30, return.model = TRUE, 
                    verbose = FALSE)
      return(obj)
    })
  }
  
  # Prepare reference data structure for MapQuery
  refdata <- setNames(as.list(idents), idents)
  
  # Perform bidirectional mapping
  mapped_objects <- list()
  
  for (direction in 1:2) {
    # Set reference and query based on direction
    # Direction 1: obj1 is reference, obj2 is query
    # Direction 2: obj2 is reference, obj1 is query
    reference <- objs[[ifelse(direction == 1, 1, 2)]]
    query <- objs[[ifelse(direction == 1, 2, 1)]]
    
    # Find anchors and transfer labels
    anchors <- FindTransferAnchors(
      reference = reference, 
      query = query, 
      reference.reduction = "pca", 
      dims = 1:30, 
      k.filter = NA
    )
    
    mapped_objects[[direction]] <- MapQuery(
      anchorset = anchors, 
      reference = reference, 
      query = query, 
      refdata = refdata, 
      reference.reduction = "pca", 
      reduction.model = "umap"
    )
  }
  
  return(mapped_objects)
}


#' MapObject
#'
#' Transfer cell type labels from a reference Seurat object to a query object using 
#' anchor-based mapping. Performs unidirectional label transfer with optional normalization.
#'
#' @param seurat_obj1 Seurat object serving as reference with known cell type labels
#' @param seurat_obj2 Seurat object serving as query to receive predicted labels
#' @param idents Character, name of metadata column in reference containing labels to 
#'   transfer (e.g., "subclass", "type", "cluster")
#' @param assay Character, assay for anchor finding: "SCT" runs SCTransform v2 (default), 
#'   "integrated" uses existing integrated assay
#' @param do.norm Logical, whether to normalize objects before mapping (default: TRUE). 
#'   Set to FALSE if objects already normalized and have PCA/UMAP reductions
#'
#' @return Query Seurat object (seurat_obj2) with added metadata:
#' \itemize{
#'   \item predicted.[idents]: Transferred labels from reference
#'   \item predicted.[idents].score: Prediction confidence scores (0-1)
#'   \item ref.umap: Query cells projected into reference UMAP coordinates
#' }
#'
#' @details
#' **Workflow:**
#' \enumerate{
#'   \item Optionally normalize both objects based on do.norm parameter
#'   \item Run PCA (30 components) and UMAP with return.model = TRUE
#'   \item Find transfer anchors between reference and query using PCA space
#'   \item Transfer labels from reference to query via MapQuery
#'   \item Project query cells into reference UMAP coordinates for visualization
#' }
#'
#' **Normalization options:**
#' - **do.norm = TRUE**: Runs SCTransform v2 or PCA+UMAP depending on assay
#'   - Use when objects contain raw/unnormalized counts
#'   - Required for first-time mapping of new objects
#'   
#' - **do.norm = FALSE**: Skips normalization, uses existing reductions
#'   - Use when objects already have SCT assay and PCA/UMAP reductions
#'   - Saves computation when mapping same objects multiple times
#'   - Required for objects that shouldn't be re-normalized
#'
#' **Assay selection:**
#' - **"SCT"**: Recommended for most cases. Runs SCTransform v2 normalization
#' - **"integrated"**: For pre-integrated datasets, runs PCA/UMAP on integrated assay
#'
#' **Label transfer mechanics:**
#' - Anchors found using first 30 PCs from reference object
#' - MapQuery handles k.weight automatically based on anchor count
#' - Confidence scores indicate prediction reliability (higher = more confident)
#' - Low scores may indicate novel cell types or poor reference coverage
#'
#' **Key differences from MapObjects:**
#' - **MapObject**: Unidirectional (returns 1 object), optional normalization
#' - **MapObjects**: Bidirectional (returns 2 objects), always normalizes, uses k.filter = NA
#'
#' @export
#' @family mapping
#'
#' @examples
#' \dontrun{
#'   # Cross-species mapping with normalization
#'   obj.opossum.mapped <- MapObject(
#'     seurat_obj1 = obj.mouse,
#'     seurat_obj2 = obj.opossum,
#'     idents = "subclass",
#'     assay = "SCT",
#'     do.norm = TRUE
#'   )
#'   
#'   # Examine prediction quality
#'   table(obj.opossum.mapped$predicted.subclass)
#'   summary(obj.opossum.mapped$predicted.subclass.score)
#'   
#'   # Visualize in reference space
#'   DimPlot(obj.opossum.mapped, reduction = "ref.umap", 
#'           group.by = "predicted.subclass", label = TRUE)
#'   
#'   # Visualize prediction confidence
#'   FeaturePlot(obj.opossum.mapped, reduction = "ref.umap",
#'               features = "predicted.subclass.score")
#'   
#'   # Within-species mapping without re-normalization
#'   # (objects already have SCT normalization and reductions)
#'   obj.test.mapped <- MapObject(
#'     seurat_obj1 = obj.train,
#'     seurat_obj2 = obj.test,
#'     idents = "cluster",
#'     assay = "SCT",
#'     do.norm = FALSE
#'   )
#'   
#'   # Generate confusion matrix
#'   confusion <- table(
#'     true = obj.test.mapped$cluster,
#'     predicted = obj.test.mapped$predicted.cluster
#'   )
#'   
#'   # Multiple label transfer
#'   obj.mapped.multi <- MapObject(
#'     seurat_obj1 = obj.reference,
#'     seurat_obj2 = obj.query,
#'     idents = c("subclass", "type"),
#'     assay = "SCT"
#'   )
#'   # Now has predicted.subclass, predicted.type, and their scores
#' }
MapObject <- function(seurat_obj1, seurat_obj2, idents, assay = "SCT", do.norm = TRUE) {
  
  library(Seurat)
  
  objs <- list(seurat_obj1, seurat_obj2)
  
  # Normalize and prepare objects if requested
  if (do.norm) {
    if (assay == "integrated") {
      objs <- lapply(objs, function(obj) {
        obj <- RunPCA(obj, npcs = 30, assay = "integrated", verbose = FALSE) %>%
               RunUMAP(reduction = "pca", dims = 1:30, assay = "integrated", 
                      return.model = TRUE, verbose = FALSE)
        return(obj)
      })
    } else if (assay == "SCT") {
      objs <- lapply(objs, function(obj) {
        obj <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = FALSE, 
                          verbose = FALSE) %>%
               RunPCA(npcs = 30, verbose = FALSE) %>%
               RunUMAP(reduction = "pca", dims = 1:30, return.model = TRUE, 
                      verbose = FALSE)
        return(obj)
      })
    }
  }
  
  # Prepare reference data structure for MapQuery
  refdata <- setNames(as.list(idents), idents)
  
  # Set reference and query
  reference <- objs[[1]]
  query <- objs[[2]]
  
  # Find anchors and transfer labels
  anchors <- FindTransferAnchors(
    reference = reference, 
    query = query, 
    reference.reduction = "pca", 
    dims = 1:30
  )
  
  obj.mapped <- MapQuery(
    anchorset = anchors, 
    reference = reference, 
    query = query, 
    refdata = refdata, 
    reference.reduction = "pca", 
    reduction.model = "umap"
  )
  
  return(obj.mapped)
}


#' PlotMapping
#'
#' Generate comprehensive visualization of label transfer results from MapObjects.
#' Creates multiple plot types showing predicted labels, confidence scores, and 
#' confusion matrices for both mapping directions.
#'
#' @param objs List of 2 mapped Seurat objects from MapObjects output
#' @param idents Character vector of metadata columns to visualize (default: c("subclass", "type"))
#' @param ident.order Optional character vector specifying order of labels in confusion matrices
#' @param title.key Character, metadata column to use for plot titles (default: "species")
#'
#' @return NULL (plots printed to graphics device). Generates the following plots for each 
#'   object and ident:
#' \itemize{
#'   \item Original UMAP with true labels
#'   \item Reference UMAP with predicted labels  
#'   \item Reference UMAP with prediction scores
#'   \item Confusion matrix (true vs predicted)
#'   \item Mapping quality heatmap (average scores per label)
#' }
#'
#' @details
#' Designed for interactive exploration of bidirectional mapping results. For each mapped
#' object, visualizes labels in both the object's native UMAP space and projected into
#' reference UMAP space. Confusion matrices show classification accuracy, while quality
#' heatmaps highlight labels with poor transfer.
#'
#' All plots use consistent axis limits (-18 to 18) and equal aspect ratios for direct
#' comparison across conditions.
#'
#' @export
#' @family mapping
#'
#' @examples
#' \dontrun{
#'   # After bidirectional mapping
#'   mapped <- MapObjects(obj.mouse, obj.opossum, idents = "subclass")
#'   
#'   # Generate all diagnostic plots
#'   PlotMapping(mapped, idents = "subclass")
#'   
#'   # Custom label order for confusion matrices
#'   PlotMapping(
#'     objs = mapped,
#'     idents = c("subclass", "type"),
#'     ident.order = c("IT_A", "IT_B", "L5PT", "L6CT")
#'   )
#' }
PlotMapping <- function(objs, idents = c("subclass", "type"), ident.order = NULL, title.key = "species") {
  
  library(Seurat)
  library(ggplot2)
  
  for (obj.idx in 1:2) {
    obj <- objs[[obj.idx]]
    
    # Get labels from the other object for confusion matrix levels
    other.idx <- setdiff(c(1, 2), obj.idx)
    
    for (ident in idents) {
      ident.levels <- as.character(unique(objs[[other.idx]][[ident, drop = TRUE]]))
      
      # Generate plots for both UMAP spaces and label types
      for (reduction in c("umap", "ref.umap")) {
        for (label.type in c("%s", "predicted.%s", "predicted.%s.score")) {
          
          label.col <- sprintf(label.type, ident)
          plot.title <- paste(obj[[title.key]][1, 1], label.col, "on", reduction)
          
          if (label.type == "predicted.%s.score") {
            # Feature plot for prediction scores
            plot <- FeaturePlot(obj, label.col, reduction = reduction) + 
                    xlim(-18, 18) + ylim(-18, 18) + coord_equal() +
                    ggtitle(plot.title)
          } else {
            # DimPlot for categorical labels
            plot <- DimPlot(obj, reduction = reduction, group.by = label.col, 
                          label = TRUE, raster = FALSE) + 
                    NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal() +
                    ggtitle(plot.title)
          }
          
          print(plot)
        }
      }
      
      # Confusion matrix
      confusion.plot <- PlotMappedLabelsHeatmap(
        obj, 
        ident, 
        ident.levels, 
        normalize = "row", 
        ident.order = ident.order
      )
      print(confusion.plot)
      
      # Mapping quality heatmap
      quality.plot <- PlotMappingQualityHeatmap(
        obj, 
        ident, 
        ident.levels, 
        sprintf("predicted.%s.score", ident), 
        ident.order = ident.order
      )
      print(quality.plot)
    }
  }
}