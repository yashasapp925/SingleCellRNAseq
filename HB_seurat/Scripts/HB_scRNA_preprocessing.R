# Load necessary libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(DoubletFinder)
library(Azimuth)

# Load data
samples <- list.dirs(path = 'GSE180665_RAW/', recursive = F, full.names = F)

for (s in samples) {
  sample_name <- gsub('_filtered_feature_bc_matrix', '', s)
  counts <- ReadMtx(mtx = paste0('GSE180665_RAW/', s, '/matrix.mtx.gz'),
                    features = paste0('GSE180665_RAW/', s, '/features.tsv.gz'),
                    cells = paste0('GSE180665_RAW/', s, '/barcodes.tsv.gz'))
  assign(sample_name, CreateSeuratObject(counts = counts))
}

# QC loop - filter low quality cells, find clusters, and remove doublets
seurat_list <- list(
  "HB17_background" = HB17_background,
  "HB17_PDX" = HB17_PDX,
  "HB17_tumor" = HB17_tumor,
  "HB30_PDX" = HB30_PDX,
  "HB30_tumor" = HB30_tumor,
  "HB53_background" = HB53_background,
  "HB53_tumor" = HB53_tumor
)

results_list <- list()

for (name in names(seurat_list)) {
  obj <- seurat_list[[name]]
  
  # Calculating mitochondrial DNA % and removing low quality cells
  obj$mtPercent <- PercentageFeatureSet(obj, pattern = '^MT-')
  obj.filt <- subset(obj, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mtPercent < 10)
  
  # Standard analysis workflow
  obj.filt <- NormalizeData(object = obj.filt)
  obj.filt <- FindVariableFeatures(object = obj.filt)
  obj.filt <- ScaleData(object = obj.filt)
  obj.filt <- RunPCA(object = obj.filt)
  
  obj.filt <- FindNeighbors(object = obj.filt, dims = 1:20)
  obj.filt <- FindClusters(object = obj.filt)
  obj.filt <- RunUMAP(object = obj.filt, dims = 1:20)
  
  # Determine cell annotations
  obj.filt <- RunAzimuth(obj.filt, reference = "HumanLiverRef/")
  
  # Determine optimal pK value for doublet detection
  sweep.res.list_hb <- paramSweep(obj.filt, PCs = 1:20, sct = FALSE)
  sweep.stats_hb <- summarizeSweep(sweep.res.list_hb, GT = FALSE)
  bcmvn_hb <- find.pK(sweep.stats_hb)
  pK <- bcmvn_hb %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  # Homotypic doublet proportion estimate
  annotations <- obj.filt@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(obj.filt@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Find and remove doublets
  obj.filt <- doubletFinder(obj.filt, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  classification <- colnames(obj.filt@meta.data)[ncol(obj.filt@meta.data)]
  obj <- obj.filt[, obj.filt@meta.data[[classification]] == "Singlet"]
  
  results_list[[name]] <- obj
}

# Merge data
HB17_background <- results_list[["HB17_background"]]
HB17_PDX <- results_list[["HB17_PDX"]]
HB17_tumor <- results_list[["HB17_tumor"]]
HB30_PDX <- results_list[["HB30_PDX"]]
HB30_tumor <- results_list[["HB30_tumor"]]
HB53_background <- results_list[["HB53_background"]]
HB53_tumor <- results_list[["HB53_tumor"]]

merged_data <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, HB53_background, HB53_tumor),
                     add.cell.ids = c("HB17_background", "HB17_PDX", "HB17_tumor", "HB30_PDX", "HB30_tumor", "HB53_background", "HB53_tumor"),
                     project = 'HB')

# Save data
SaveSeuratRds(merged_data, file = 'HB_merged_samples.rds')