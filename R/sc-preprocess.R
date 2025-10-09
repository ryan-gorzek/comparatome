#' PreprocessData
#'
#' Preprocess multiple 10X-like directories containing sn or scRNA-seq data.
#' Each directory should exist within data_path and its name be included in sample_IDs.
#' Directories must meet the criteria for Seurat v4's Read10X.
#' @param data_path Path to the directory containing sample_IDs folder with RNA-seq data in 10X format.
#' @param sample_IDs A character vector of folder names within data_path, each containing files meeting the requirements for Seurat v4's Read10X. Encoded as 'sample' in the metadata of the resulting Seurat object.
#' @param project_name String passsed to the 'project' argument of CreateSeuratObject.
#' @param mapping_path Path to a TXT file containing cross-species gene orthologs. See documentation for MapGenes.
#' @param gene.column Column number passed to the 'gene.column' argument of Read10X.
#' @return List of pre- and post-filtered (standard cell and gene count criteria) objects and gene lists from all specified sample_IDs, with sample and project identifiiers, and Scrublet metadata.
#' @export
#' @family preprocess
#' @examples
#' \dontrun{
#'  data_path <- paste0(dir.list$central, "samples/")
#'  sample_IDs <- c('OpossumV1-3A', 'OpossumV1-3B', 'OpossumV1-4A', 'OpossumV1-4B')
#'  mapping_path <- paste0(dir.list$central, "Opossum_Mouse_GeneMapping_EnsemblBioMart.txt")
#'
#'  data <- PreprocessData(data_path, sample_IDs, "Opossum_V1", mapping_path)
#'
#'  obj.opossum <- data$obj
#' }
PreprocessData <- function(data_path, sample_IDs, project_name, mapping_path = NA, gene.column = 2) {
  
  # Load data for each sample
  print("Loading data...")
  objs <- c()
  for (sample in sample_IDs) {
    
    temp.obj.path <- paste0(data_path, sample)
    temp.obj.data <- Read10X(temp.obj.path, gene.column = gene.column)
    temp.obj <- CreateSeuratObject(counts = temp.obj.data, project = project_name)
    temp.obj$sample <- sample
    temp.obj <- scrublet_R(seurat_obj = temp.obj)
    objs <- append(objs, temp.obj)
    
  }
  # Combine samples into single Seurat object
  if (length(objs) > 1) {
    obj <- merge(objs[[1]], y = objs[2:length(objs)], add.cell.ids = sample_IDs, project = project_name)
  } else { 
    obj <- objs[[1]] 
  }
  # Store original gene names
  genes <- list()
  genes$pre.map <- rownames(obj)
  # Map gene names/IDs if mapping_path is specified
  if !is.na(mapping_path) {
    obj <- MapGenes(obj, mapping_path)
    genes$post.map <- rownames(obj)
  }
  # Filter cells and genes with standard critera (based on Cheng et al., Cell 2022)
  print("Filtering cells and genes...")
  obj.prefilt <- obj
  cell_mask <- Reduce(intersect,list(WhichCells(obj, expression = nFeature_RNA > 700),
                                     WhichCells(obj, expression = nFeature_RNA < 6500),
                                     WhichCells(obj, expression = nCount_RNA < 40000)))
  gene_mask <- rownames(obj)[Matrix::rowSums(obj[["RNA"]]@counts > 0) > 8]
  obj <- subset(obj, features = gene_mask, cells = cell_mask)
  # Store post-filtered gene names
  genes$post.filt <- rownames(obj)
  # Return list containing various objects
  return(list(obj = obj, obj.prefilt = obj.prefilt, genes = genes))
}

#' MapGenes
#'
#' Map gene names between species given a Seurat object to be mapped, and an orthology table.
#' Orthology tables are most easily obtained from Ensembl BioMart. For an example, see:
#' github.com/ryan-gorzek/opossum-V1-omics/blob/main/central/Opossum_Mouse_GeneMapping_EnsemblBioMart.txt
#' Columns must be: 
#'    1. Gene name (mapped from)
#'    2. Gene ID (matching column 1)
#'    3. Gene name (mapped to)
#'    4. Gene ID (matching column 3)
#' @param obj Seurat object with gene names/IDs to be mapped.
#' @param mapping_path Path to a file containing a gene orthology table. See MapGenes description for formatting details.
#' @param use_ids T/F parameter that specifies whether to rely on gene IDs for matching. Set to FALSE if you've built a custom mapping file without gene IDs.
#' @return Seurat object with mapped gene names/IDs.
#' @export
#' @family sc-preprocess
#' @examples
#' \dontrun{
#'  seurat.obj.opossum <- MapGenes(seurat.obj.opossum, ../Opossum_Mouse_GeneMapping_EnsemblBioMart.txt)
#' }
MapGenes <- function(obj, mapping_path, use_ids = TRUE) {
  
  # Read the orthology table
  genes.mapping <- read.csv(mapping_path)
  # Remove paralogs (keeping the first)
  para.idx <- genes.mapping[, 1] %in% unique(genes.mapping[duplicated(genes.mapping[, 1]), 1])
  genes.mapping <- genes.mapping[!para.idx,]
  # Store columns of interest separately
  genes.mapping.self <- as.list(genes.mapping[, 1])
  ids.mapping.self <- as.list(genes.mapping[, 2])
  genes.mapping.other <- as.list(genes.mapping[, 3])
  genes.self <- rownames(obj)
  # Loop over genes (from 'other', which is mapped to) that have orthologs (i.e., exist in table)
  for (gene in genes.mapping.other) {
    # Find all entries that match current gene
    idx.other <- which(genes.mapping.other %in% gene)
    # Proceed only if there is a 1:1 mapping between 'other' gene and 'self' gene
    if (length(idx.other) == 1) {
      gene.self <- genes.mapping.self[idx.other]
      id.self <- ids.mapping.self[idx.other]
      # If gene is not named in 'self' (e.g., novel gene only known by ID) or use_IDs is T, map based on ID
      if ((gene.self == "") | (use_ids == TRUE)) {
        idx.self <- which(genes.self %in% id.self)
        genes.self[idx.self] <- gene
        # Otherwise, use gene name if it exists
      } else if (gene.self != "") {
        idx.self <- which(genes.self %in% gene.self)
        genes.self[idx.self] <- gene
      }
    }
  }
  # Rebuild Seurat object
  print("Rebuilding object...")
  obj.df <- as.data.frame(as.matrix(obj[["RNA"]]@counts))
  rownames(obj.df) <- genes.self
  obj.temp <- CreateSeuratObject(counts = obj.df, meta.data = obj[[]])
  obj <- obj.temp
  return(obj)
}


#' NormalizeAndPCA
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the preprocess family.
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
#' Part of the preprocess family.
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
