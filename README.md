# ABCT
ABCT(Anchor-based Cell Typer)

![Graphical Abstract](images/graphical_abstract.png)


## Overview

ABCT is a hybrid cell type annotation method designed for spatial omics data, combining the simplicity of marker-based annotation with the reliability of profile-based methods. It calculates cell type scores using a list of marker genes for each cell type and incorporates spatial information by considering the expression of neighboring cells. This enables the identification of anchor cells that represent each cell type, creating a profile that is then used to annotate the entire dataset.

The key advantages of ABCT include:

- No complex integration with scRNA-seq datasets required
- Support for novel cell types defined by the user
- Effective utilization of spatial information for accurate cell type classification
- Clear differentiation between malignant and normal epithelial cells in tumor samples

ABCT is compatible with a wide range of spatial technologies (e.g., 10x Xenium, MERFISH, CosMX, CODEX) and scales efficiently for large datasets. For more information, check out:

- the [paper](https://www.nature.com/articles/s41588-024-01664-3), (논문 링크)
- a tutorial on ABCT, a set of [vignettes](https://prabhakarlab.github.io/Banksy) showing basic usage, usage compatibility with Seurat ([here](https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md) and [here](https://satijalab.org/seurat/articles/visiumhd_analysis_vignette#identifying-spatially-defined-tissue-domains)),  (tutorial)
- a [Zenodo archive](https://zenodo.org/records/10258795) containing scripts to reproduce the analyses in the paper (데이터 다운로드 페이지?)


## Quick start
### Loading Required Libraries
Ensure that all necessary libraries are installed and loaded.

```{r libraries}
library(Seurat)
library(SeuratWrappers)
library(Banksy)
library(harmony)
library(dplyr)
library(scales)
library(EnvStats)
library(stringr)
library(matrixStats)
library(InSituType)
library(UCell)
library(ggplot2)
library(MASS)
library(pracma)
source("./ABCT_final.r")
```
