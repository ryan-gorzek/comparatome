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
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
#' @param folder_path (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
make_folder <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    # Create the folder
    dir.create(folder_path)
  }
}


#' strip_if_contains
#'
#' Auto-generated roxygen skeleton for comparatome.
#' Part of the helpers family.
#' @param vector (auto) parameter
#' @param search (auto) parameter
#' @param strip (auto) parameter
#' @return (auto) value; see function body.
#' @keywords internal
#' @family helpers
#' @examples
#' \dontrun{
#'  # Example usage will be added
#' }
strip_if_contains <- function(vector, search, strip) {
  # Initialize an empty vector to store results
  result_vector <- vector
  
  # Loop through each element of the vector
  for (i in seq_along(vector)) {
    if (grepl(search, vector[i])) {
      result_vector[i] <- paste0("U#", sub(paste0("^", strip), "", vector[i]))
    }
  }
  
  return(result_vector)
}
