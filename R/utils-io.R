#' URLencodeNCBI
#'
#' Encode a string (filename) with NCBI URL encoding for programmatic download from Gene Expression Omnibus.
#' R's standard URLencode function does not suffice.
#' @param x String (typically filename) for programmatic download from NCBI Gene Expression Omnibus.
#' @return String encoded for NCBI URLs.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", samples.opossum[[sm]],
#'                "&format=file&file=", URLencodeNCBI(fname))
#'  system(paste0("wget -O ", fpath, " ", url))
#' }
URLencodeNCBI <- function(x) {
  b <- charToRaw(enc2utf8(x))
  out <- vapply(b, function(xx) {
    if ((xx >= 0x30 && xx <= 0x39) || (xx >= 0x41 && xx <= 0x5A) || (xx >= 0x61 && xx <= 0x7A)) {
      rawToChar(as.raw(xx))
    } else {
      sprintf("%%%02X", as.integer(xx))
    }
  }, character(1))
  paste0(out, collapse = "")
}


#' SavePNGandSVG
#'
#' Saves a ggplot object as both PNG and SVG files to a specified directory.
#' Automates exporting a plot in two formatsusing the same filename and output path. 
#' Loops through the formats and calls [ggplot2::ggsave()] for each.
#'
#' @param p A ggplot object to be saved.
#' @param fpath A string specifying the directory path where files will be saved. Must include a trailing slash (e.g., `"figures/"`).
#' @param fname A string specifying the base filename (without file extension).
#' @return Invisibly returns `NULL`. Files are written to disk as a side effect.
#' @keywords internal
#' @family helpers
#' @importFrom ggplot2 ggsave
#' @examples
#' \dontrun{
#'  obj.opossum <- ClusterWithSCT(obj.opossum, 1)
#'  clust.plots <- PlotClusters(obj.opossum)
#'  save.plots <- list("dimplot_cluster" = "S1E_All-Cluster-DimPlot",
#'                     "dimplot_sample" = "S1F_All-Sample-DimPlot",
#'                     "barplot_sample" = "S1G_All-Sample-BarPlot",
#'                     "barplot_doublet" = "S1H_All-Doublet-BarPlot")
#'  for (sp in names(save.plots)) {
#'    SavePNGandSVG(clust.plots[[sp]], dir.list$figS1$plots, save.plots[[sp]])
#'  }
#' }
SavePNGandSVG <- function(p, fpath, fname) {
  for (ft in c(".png", ".svg")) {
    ggsave(paste0(fpath, fname, ft), plot = p)
  }
}


#' make_folder
#'
#' Create a directory if it does not already exist.
#' Convenience wrapper around dir.exists() and dir.create() for safe directory creation.
#'
#' @param folder_path String specifying the path to the directory to create.
#'   Can be a relative or absolute path.
#'
#' @return NULL (invisibly). Creates directory as a side effect.
#'
#' @details
#' This function checks if the specified directory exists. If it does not exist,
#' the directory is created. If it already exists, no action is taken and no
#' warning is generated.
#'
#' Useful for ensuring output directories exist before saving files, preventing
#' errors from attempting to write to non-existent directories.
#'
#' @keywords internal
#' @family helpers
#'
#' @examples
#' \dontrun{
#'  # Ensure plot directory exists before saving
#'  make_folder("figures/dimplots/")
#'  ggsave("figures/dimplots/cluster_umap.png", plot = p)
#'  
#'  # Create nested directory structure
#'  make_folder("output/markers/IT_A/")
#' }
make_folder <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    # Create the folder
    dir.create(folder_path, recursive = TRUE)
  }
}


#' strip_if_contains
#'
#' Strip a specific pattern from strings that contain a search term.
#' Used to clean gene names by removing Ensembl ID prefixes when present.
#'
#' @param vector Character vector of strings to process (e.g., gene names).
#' @param search String pattern to search for within each element. If this pattern
#'   is found, the strip operation is performed.
#' @param strip String pattern to remove from the beginning of matching elements.
#'   Typically a longer, more specific version of the search pattern.
#'
#' @return Character vector with the same length as input. Elements containing the
#'   search pattern have the strip pattern removed from the beginning and are prefixed
#'   with "U#". Elements not containing the search pattern are returned unchanged.
#'
#' @details
#' The function processes each element of the input vector:
#' 1. Checks if the element contains the search pattern
#' 2. If yes, removes the strip pattern from the beginning using gsub with anchored regex
#' 3. Prefixes the result with "U#" to indicate an unnamed/unknown gene
#' 4. If no, returns the element unchanged
#'
#' This is particularly useful for cleaning Ensembl gene identifiers in plots.
#' For example, "ENSMUSG00000012345" might be stripped to "U#12345" for cleaner
#' visualization while preserving the ability to identify the gene.
#'
#' @keywords internal
#' @family helpers
#'
#' @examples
#' \dontrun{
#'  # Clean mouse Ensembl IDs
#'  genes <- c("Slc17a6", "ENSMUSG00000025400", "Gad1", "ENSMUSG00000070880")
#'  strip_if_contains(genes, "ENSMUSG", "ENSMUSG000000")
#'  # Returns: c("Slc17a6", "U#25400", "Gad1", "U#70880")
#'  
#'  # Clean opossum Ensembl IDs
#'  genes <- c("SLC17A6", "ENSMODG00000012345", "GAD1")
#'  strip_if_contains(genes, "ENSMODG", "ENSMODG000000")
#'  # Returns: c("SLC17A6", "U#12345", "GAD1")
#' }
strip_if_contains <- function(vector, search, strip) {
  # Initialize result vector as copy of input
  result_vector <- vector
  
  # Loop through each element of the vector
  for (i in seq_along(vector)) {
    if (grepl(search, vector[i])) {
      # If search pattern found, remove strip pattern and add U# prefix
      result_vector[i] <- paste0("U#", sub(paste0("^", strip), "", vector[i]))
    }
  }
  
  return(result_vector)
}
