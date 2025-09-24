#' comparatome: Comparative single-cell and cross-species transcriptomics utilities
#'
#' Tools built on Seurat for preprocessing, clustering, mapping across datasets/species,
#' cross-species DE comparisons, and diagnostic visualization.
#'
#' @docType package
#' @name comparatome
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP FindNeighbors FindClusters DimPlot FeaturePlot VlnPlot DoHeatmap FindTransferAnchors IntegrateData TransferData AverageExpression VariableFeatures
#' @importFrom SeuratObject GetAssayData DefaultAssay Idents
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_histogram geom_bar geom_tile scale_fill_manual scale_color_manual theme theme_minimal labs xlab ylab coord_fixed geom_text
#' @importFrom dplyr mutate select filter arrange group_by summarise left_join inner_join bind_rows bind_cols distinct rename
#' @importFrom tidyr pivot_longer pivot_wider drop_na
#' @importFrom Matrix rowSums colSums sparseMatrix
#' @importFrom matrixStats rowVars colVars
#' @importFrom stringr str_detect str_replace str_remove str_to_lower
#' @importFrom gridExtra grid.arrange
#' @importFrom stats cor dist prcomp setNames median sd lm p.adjust
NULL