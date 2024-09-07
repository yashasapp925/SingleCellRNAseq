# Load necessary libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Load data
obj <- LoadSeuratRds("HB_integrated_samples.rds")

# Clustering integrated data
ElbowPlot(obj)

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:20)
obj <- FindClusters(obj, cluster.name = "cca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")

# Visualizations that demonstrate integration of data
patient_intg <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = 'Patient_ID',
)
tissue_intg <- DimPlot(
  obj, 
  reduction = 'umap.cca', 
  group.by = 'Tissue Type', 
  cols = c('red', 'green', 'blue')
)
png("Plots/integrated_data_categories.png", width = 800, height = 600)
grid.arrange(patient_intg, tissue_intg, ncol = 2, nrow = 2)
dev.off()

cca_clusters <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = 'cca_clusters',
)
png("Plots/integrated_data_clusters.png", width = 800, height = 600)
cca_clusters
dev.off()

celltype <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c('predicted.celltype.l2')
)
png("Plots/integrated_data_celltypes.png", width = 800, height = 600)
celltype
dev.off()

# Save data
SaveSeuratRds(obj, file = 'HB_integrated_clustered.rds')
