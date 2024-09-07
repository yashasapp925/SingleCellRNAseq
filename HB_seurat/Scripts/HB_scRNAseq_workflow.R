# Load necessary libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Load data
merged_data <- LoadSeuratRds("HB_merged_samples.rds")

# Edit metadata
merged_data$sample <- rownames(merged_data@meta.data)
merged_data@meta.data <- separate(merged_data@meta.data, col = 'sample', into = c('Patient_ID', 'Tissue Type', 'Barcode'), sep = '_')

# Standard analysis without integration
merged_data <- NormalizeData(object = merged_data)
merged_data <- FindVariableFeatures(object = merged_data)
merged_data <- ScaleData(object = merged_data)
merged_data <- RunPCA(object = merged_data)

ElbowPlot(merged_data)

merged_data <- FindNeighbors(object = merged_data, dims = 1:20)
merged_data <- FindClusters(object = merged_data)
merged_data <- RunUMAP(object = merged_data, dims = 1:20)

# Plots
patient_plt <- DimPlot(merged_data, reduction = 'umap', group.by = 'Patient_ID')
tissue_plt <- DimPlot(merged_data, reduction = 'umap', group.by = 'Tissue Type', cols = c('red', 'green', 'blue'))
png("Plots/unintegrated_data.png", width = 800, height = 600)
grid.arrange(patient_plt, tissue_plt, ncol = 2, nrow = 2)
dev.off()

# Integration
obj <- IntegrateLayers(
  object = merged_data, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

# Save data
SaveSeuratRds(obj, file = 'HB_integrated_samples.rds')
