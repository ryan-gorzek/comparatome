#' PreprocessData
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the preprocess family.
#' @param sample_IDs (auto) parameter
#' @param data_path (auto) parameter
#' @param project_name (auto) parameter
#' @param mapping_path (auto) parameter
#' @param gene.column (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family preprocess
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PreprocessData <- function(sample_IDs, data_path, project_name, mapping_path, gene.column = 2) {
  
  # Load the data.
  print("Loading 10x data...")
  objs <- c()
  genes <- list()
  
  for (sample in sample_IDs) {
    
    temp.obj.path <- paste(data_path, sample, "/outs/filtered_feature_bc_matrix/", sep = "")
    temp.obj.data <- Read10X(temp.obj.path, gene.column = gene.column)
    temp.obj <- CreateSeuratObject(counts = temp.obj.data, project = project_name)
    temp.obj$sample <- sample
    temp.obj <- scrublet_R(seurat_obj = temp.obj)
    objs <- append(objs, temp.obj)
    
  }
  
  if (length(objs) > 1) {
  obj <- merge(objs[[1]], y = objs[2:length(objs)], add.cell.ids = sample_IDs, project = project_name)
  } else { obj <- objs[[1]] }
  genes$pre.map <- rownames(obj)
  
  # Map gene names/IDs if specified.
  print("Mapping genes...")
  if (!is.na(mapping_path)) {
    genes.mapping <- read.csv(mapping_path)
    para.idx <- genes.mapping$Gene.stable.ID %in% unique(genes.mapping$Gene.stable.ID[duplicated(genes.mapping$Gene.stable.ID)])
    genes.mapping <- genes.mapping[!para.idx,]
    genes.mapping.other <- as.list(genes.mapping[, 3])
    genes.mapping.self <- as.list(genes.mapping[, 1])
    ids.mapping.self <- as.list(genes.mapping[, 2])
    genes.self <- rownames(obj)
    for (gene in genes.mapping.other) {
      
      idx.other <- which(genes.mapping.other %in% gene)
      
      if (length(idx.other) == 1) {
        
        gene.self <- genes.mapping.self[idx.other]
        id.self <- ids.mapping.self[idx.other]
        
        if (gene.self == "") {
          
          idx.self <- which(genes.self %in% id.self)
          genes.self[idx.self] <- gene
          
        } else {
          
          idx.self <- which(genes.self %in% gene.self)
          genes.self[idx.self] <- gene
          
        }
      }
    }
  }
  else { genes.self <- rownames(obj) }
  
  # Rebuild Seurat object.
  print("Rebuilding object...")
  obj.df <- as.data.frame(as.matrix(obj[["RNA"]]@counts))
  rownames(obj.df) <- genes.self
  obj.temp <- CreateSeuratObject(counts = obj.df, meta.data = obj[[]])
  obj <- obj.temp
  genes$post.map <- rownames(obj)
  
  # Filter cells and genes.
  print("Filtering cells and genes...")
  obj.prefilt <- obj
  cell_mask <- Reduce(intersect,list(WhichCells(obj, expression = nFeature_RNA > 700),
                                     WhichCells(obj, expression = nFeature_RNA < 6500),
                                     WhichCells(obj, expression = nCount_RNA < 40000)))
  
  gene_mask <- rownames(obj)[Matrix::rowSums(obj[["RNA"]]@counts > 0) > 8]
  
  obj <- subset(obj, features = gene_mask, cells = cell_mask)
  
  genes$post.filt <- rownames(obj)
  
  return(list(obj = obj, obj.prefilt = obj.prefilt, genes = genes))
  
}


#' NormalizeAndPCA
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the preprocess family.
#' @param obj (auto) parameter
#' @param nfeatures (auto) parameter
#' @param npcs (auto) parameter
#' @param features (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family preprocess
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
NormalizeAndPCA <- function(obj, nfeatures = 3000, npcs = 30, features = NA) {
  
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  if (all(is.na(features))) {
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  } else {
    VariableFeatures(obj) <- features
  }
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = npcs)
  return(obj)
  
}


#' PlotQC
#'
#' Auto-generated roxygen skeleton for comparatome.
Part of the preprocess family.
#' @param data (auto) parameter
#' @return (auto) value; see function body.
#' @export
#' @family preprocess
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
PlotQC <- function(data) {
  
  # Create a custom theme for the violin plots
  custom_theme_vln <- theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = "none", # Remove legends
    axis.title.x = element_blank(), # Remove x-axis labels
    axis.text.x = element_text(size = 8) # Remove x-axis tick labels
  )
  
  # Create a custom theme for the scatter plots
  custom_theme_scatter <- theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.position = "none", # Remove legends
    axis.title.x = element_text(size = 8), # Add x-axis labels back
    axis.text.x = element_text(size = 6) # Add x-axis tick labels back
  )
  
  # Create individual plots for pre-filtered data with the custom theme
  vln_plot1_pre <- VlnPlot(data$obj.prefilt, features = c("nFeature_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 10000) + custom_theme_vln
  vln_plot2_pre <- VlnPlot(data$obj.prefilt, features = c("nCount_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 40000) + custom_theme_vln
  scatter_plot_pre <- FeatureScatter(data$obj.prefilt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE) + xlim(0, 40000) + ylim(0, 10000) + custom_theme_scatter + coord_fixed(ratio = 4)
  
  # Create individual plots for post-filtered data with the custom theme
  vln_plot1_post <- VlnPlot(data$obj, features = c("nFeature_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 10000) + custom_theme_vln
  vln_plot2_post <- VlnPlot(data$obj, features = c("nCount_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 40000) + custom_theme_vln
  scatter_plot_post <- FeatureScatter(data$obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE) + xlim(0, 40000) + ylim(0, 10000) + custom_theme_scatter + coord_fixed(ratio = 4)
  
  final_layout_feature <- vln_plot1_pre | vln_plot1_post
  final_layout_count <- vln_plot2_pre | vln_plot2_post
  final_layout_scatter <- scatter_plot_pre | scatter_plot_post
  
  # Display the final layout
  print(final_layout_feature)
  print(final_layout_count)
  print(final_layout_scatter)
  
}
