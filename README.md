# Single Cell RNAseq Projects
This repository houses projects that analyze single cell RNA-seq data using the tools Seurat and Scanpy. Different datasets are also analyzed.

## Immune Cell Project
This project involves the doublet removal, preprocessing, integration, clustering, and manual cell type labeling of several immune cell samples from atherosclerosis plaques. It uses scanpy for single-cell analysis and scvi-tools to remove doublets and integrate the samples.

### Dataset
The data for this project was acquired from the Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245373). The original paper that generated and studied this data was "Cross-species single-cell RNA sequencing reveals divergent phenotypes and activation states of adaptive immunity in human carotid and experimental murine atherosclerosis" by Horstmann et al. (10.1093/cvr/cvae154).

### Doublet removal
Using the scvi-tools SOLO model, doublets were predicted and removed from the sample datasets. SOLO is a generative model of scRNA-seq count data that is used to detect doublets, and this is done by first training a VAE model on a single sample, using that model to simulate doublets, then training a classifier to finally predict whether each cell in the sample is a doublet or not (https://docs.scvi-tools.org/en/stable/user_guide/models/solo.html).

### Cell type annotation
The annotation for this project was done manually by examining differentially expressed genes between clusters and using PangloaDB to find common immune cell type markers, which were used for inspecting the clusters.

## Hepatoblastoma Project
This project involves the integration of multiple liver cell samples from Hepatoblastoma patients. It includes a full scRNA analysis workflow with quality control steps such as removing cells with high mitochondrial DNA percentage and doublet removal. It finishes with cell clustering and an application of an automatic cell type annotation tool called Azimuth.

### Dataset
The data for this project was acquired from the Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180665).

### Cell type annotation
An important component of this project was automatic cell type anotation using the tool Azimuth (https://github.com/satijalab/azimuth). In order to annotate the cells, a liver cell reference dataset was required (https://zenodo.org/records/7770308).
