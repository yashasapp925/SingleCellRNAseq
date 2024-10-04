# Single Cell RNAseq Projects
This repository houses projects that analyze single cell RNA-seq data using the tools Seurat and Scanpy. Different datasets are also analyzed.

## Hepatoblastoma Project
This project involves the integration of multiple liver cell samples from Hepatoblastoma patients. It includes a full scRNA analysis workflow with quality control steps such as removing cells with high mitochondrial DNA percentage and doublet removal. It finishes with cell clustering and an application of an automatic cell type annotation tool called Azimuth.

### Dataset
The data for this project was acquired from the Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180665).

### Cell type annotation
An important component of this project was automatic cell type anotation using the tool Azimuth (https://github.com/satijalab/azimuth). In order to annotate the cells, a liver cell reference dataset was required (https://zenodo.org/records/7770308).
