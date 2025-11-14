# comparatome

Comparative (cross-dataset/species) single-cell and spatial transcriptomics tools built on Seurat for preprocessing, clustering, integration/mapping, gene expression comparisons, and visualization. Under active development.

## Installation

```r
# Install from GitHub
devtools::install_github("ryan-gorzek/comparatome")
```

## Features

### Single-cell Preprocessing (`sc-*`)
- **Data Loading & QC**: Load 10X data, filter cells/genes, integrate Scrublet doublet detection
- **Gene Mapping**: Cross-species ortholog mapping via Ensembl BioMart tables
- **Normalization**: SCTransform v2, log-normalization, PCA, UMAP
- **Clustering**: Leiden clustering across multiple resolutions
- **Marker Discovery**: Differential expression analysis at cluster and subclass levels
- **Visualization**: Dimensionality reduction plots, feature plots, dot plots, heatmaps

### Cross-species Analysis (`cross-*`)
- **Integration**: Integrate datasets/species using Seurat v4 CCA anchors
- **Mapping**: Transfer cell type labels between datasets/species
- **Mapping Metrics**: Quantitative assessment of cross-species correspondence
- **Confusion Matrices**: Cell type alignment quantification
- **DE Comparisons**: Cross-species differential expression curve analysis
- **DE Overlap**: Intersection analysis of marker genes across species
- **Component Analysis**: CCA feature loading visualization, elbow plot comparisons
- **Heatmaps**: Cross-species cluster overlap matrices

### Utilities (`utils-*`)
- **I/O**: PNG/SVG export, NCBI URL encoding for GEO downloads
- **Gene Lists**: Gene set operations
- **Color Management**: Consistent color schemes
- **Matrix Operations**: Sparse matrix utilities
- **Ordering**: Cluster/subclass sorting by reference hierarchies

## Function Families

Functions are organized into families for documentation:

- `preprocess`: Data loading, filtering, normalization, gene mapping
- `cluster`: Clustering and cluster visualization
- `markers`: Differential expression and marker gene identification
- `visualization`: Plotting utilities
- `mapping`: Cross-dataset/species mapping
- `integration`: Dataset integration
- `mapping-metrics`: Quantitative mapping assessment
- `confusion`: Confusion matrix generation
- `de-curves`: Differential expression comparison curves
- `de-overlap`: Marker gene overlap analysis
- `components`: Dimensionality reduction analysis
- `heatmaps`: Cross-species overlap heatmaps
- `subsampling`: Expression shuffling and subsampling
- `helpers`: I/O, color, ordering utilities

## Example Workflow

```r
# Load and preprocess data
data <- PreprocessData(
  data_path = "path/to/samples/",
  sample_IDs = c("Sample1", "Sample2"),
  project_name = "MyProject",
  mapping_path = "path/to/ortholog_table.txt"
)
obj <- data$obj

# Cluster with SCTransform
obj <- ClusterWithSCT(obj, resolutions = c(0.2, 0.5, 1.0))

# Generate cluster plots
plots <- PlotClusters(obj, group.id = "SCT_snn_res.0.5")

# Find markers
markers <- SubclassMarkerDict(
  obj,
  subclass.col = "subclass",
  save.path = "markers.rds"
)

# Cross-species mapping
obj.mapped <- MapObject(
  seurat_obj1 = obj.species1,
  seurat_obj2 = obj.species2,
  idents = c("subclass", "cluster")
)

# Visualize mapping with confusion matrix
p <- PlotMappedLabelsHeatmap(
  data = obj.mapped,
  column_name = "subclass",
  column_levels = c("IT", "Pvalb", "Sst", "Vip"),
  normalize = "row"
)
```

## Documentation

All exported functions include roxygen2 documentation. Access help:

```r
?PreprocessData
?ClusterWithSCT
?MapObject
```

## Roadmap

### In Development
- Expanded marker gene visualization tools
- Enhanced cross-species integration workflows
- Additional mapping metrics

### Future Plans

Additional branches:
- **Seurat v4**: Maintain compatibility with Seurat v4 for legacy analyses
- **Seurat v5**: Full integration with Seurat v5 features and BPCells backend
- **Scanpy**: Python implementations leveraging Scanpy, AnnData, and scikit-learn ecosystem
- **Additional Python Tools**: rapids-singlecell and other Python-based tools for large-scale analysis

These branches will enable users to select the optimal computational environment for their workflows as datasets grow to millions of cells.

## License

MIT. See `LICENSE`.
