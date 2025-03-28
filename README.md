# ABCT
ABCT(Anchor-based Cell Typer)

![Graphical Abstract](images/0_graphical_abstract.png)


## Overview

ABCT is a hybrid cell type annotation method designed for spatial omics data, combining the simplicity of marker-based annotation with the reliability of profile-based methods. It calculates cell type scores using a list of marker genes for each cell type and incorporates spatial information by considering the expression of neighboring cells. This enables the identification of anchor cells that represent each cell type, creating a profile that is then used to annotate the entire dataset.

The key advantages of ABCT include:

- No complex integration with scRNA-seq datasets required
- Support for novel cell types defined by the user
- Effective utilization of spatial information for accurate cell type classification
- Clear differentiation between malignant and normal epithelial cells in tumor samples

ABCT is compatible with a wide range of spatial technologies (e.g., 10x Xenium, MERFISH, CosMX, CODEX) and scales efficiently for large datasets. For more information, check out:

- the [paper]()
- a tutorial on ABCT, a set of [vignettes]()
- [HISSTA](https://kbds.re.kr/hissta/) containing data used in the analyses presented in the paper





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
source("./ABCT.r")
```

### Loading and Preparing Data
Load the data for analysis. Here, we demonstrate using a Seurat object. The processed data can be downloaded from [HISSTA](https://kbds.re.kr/hissta/datasetinfo?sampleIdx=19).
```{r load-data}
obj <- readRDS(paste0("/path/to/your/data/", "obj_final.rds"))
```

### Define Cell Type List and Marker Genes
We define a list of cell types and specify marker genes for the analysis.
```{r cell-types}
celltype_list <- c("Malignant", "Epithelial", "CD4T", "CD8T", "NK", "B", "Plasma",
                   "Macrophage", "Monocyte", "DC", "Mast", "Neutrophil", 
                   "Endothelial", "Fibroblast", "Unknown")
names(celltype_list) <- c(
  "#bf4040", "#FB8072", "#8DD3C7", "#FFFFB3", "#BEBADA", "#3288bd", 
  "#FDB462", "#B3DE69", "#FCCDE5", "#84cdee", "#BC80BD", "#CCEBC5", 
  "#FFED6F", "#1b9e77", "gray"
)
malignant_marker_list <- data.frame(
  cluster = "Malignant",
  gene = c("SOX9", "FGFR1", "KRAS", "MYC", "EGFR")
)
ABCT_marker_list <- read.csv("/path/to/marker_list.csv")
ABCT_marker_list$cluster <- factor(ABCT_marker_list$cluster, levels = celltype_list)
ABCT_marker_list <- ABCT_marker_list %>% arrange(cluster)
```

### Identifying Malignant Cells
Use the `FindMalignantCells` function to identify malignant cells in the data.
```{r find-malignant}
malignant_path <- "/path/to/output/malignant/"
obj <- FindMalignantCells(
  obj,
  assay = "SCT",
  ctrl_assay = "negprobes",
  marker_list = malignant_marker_list,
  use_spatial = TRUE,
  M = 1,
  lambda = 0.2,
  area = NULL,
  w_neg = 0,
  dimx = "x_global_px",
  dimy = "y_global_px",
  smooth_reduction = "spatial_pca",
  path = malignant_path
)
```
<p align="center">
  <img src="images/1_malignant_anchor.png" width="30%" />
  <img src="images/1_malignant_dim.png" width="30%" />
  <img src="images/1_malignant_image.png" width="30%" />
</p>



### Performing ABCT Classification
Subsequently, classify the non-malignant cells using the RunABCT function. By default, if a BANKSY object already exists in the Seurat object, RunABCT will use the existing BANKSY object rather than running it again. If you wish to rerun BANKSY on a subsetted object, you should remove the existing BANKSY assay before executing the function.
```{r run-abct}
subobj <- subset(obj, subset = malignant_result != "Malignant")
abct_path <- "/path/to/output/abct/"
subobj <- RunABCT(
  subobj,
  assay = "SCT",
  ctrl_assay = "negprobes",
  marker_list = ABCT_marker_list,
  color_list = names(celltype_list[celltype_list %in% ABCT_marker_list$cluster]),
  method = "quantile",
  use_spatial = TRUE,
  M = 1,
  lambda = 0.2,
  dimx = "x_global_px",
  dimy = "y_global_px",
  smooth_reduction = "spatial_pca",
  path = abct_path
)
```
<p align="center">
  <img src="images/2_abct_anchor.png" width="30%" />
  <img src="images/2_abct_dim.png" width="30%" />
  <img src="images/2_abct_image.png" width="30%" />
</p>



### Updating Metadata
Finally, update the metadata in the original object with the ABCT classification results.
```{r update-metadata}
obj <- update_metadata(obj, subobj, celltype_list)
```

<br>
<br>
<details>
<summary><strong>Session Info</strong></summary>
  
```r
sessionInfo()
#> R version 4.3.2 (2023-10-31)
#> Platform: x86_64-conda-linux-gnu (64-bit)
#> Running under: Rocky Linux 8.10 (Green Obsidian)
#> 
#> Matrix products: default
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=ko_KR.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=ko_KR.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=ko_KR.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=ko_KR.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Asia/Seoul
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] pracma_2.4.4         MASS_7.3-60          ggplot2_3.5.1       
#>  [4] UCell_2.6.2          InSituType_2.0       matrixStats_1.3.0   
#>  [7] stringr_1.5.1        EnvStats_2.8.1       scales_1.3.0        
#> [10] dplyr_1.1.4          harmony_1.2.0        Rcpp_1.0.12         
#> [13] Banksy_0.99.13       SeuratWrappers_0.3.5 Seurat_5.0.1        
#> [16] SeuratObject_5.0.1   sp_2.1-4            
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppHungarian_0.3           RcppAnnoy_0.0.22           
#>   [3] splines_4.3.2               later_1.3.2                
#>   [5] bitops_1.0-7                tibble_3.2.1               
#>   [7] R.oo_1.26.0                 polyclip_1.10-6            
#>   [9] fastDummies_1.7.3           lifecycle_1.0.4            
#>  [11] aricode_1.0.3               globals_0.16.3             
#>  [13] lattice_0.22-6              SnowballC_0.7.1            
#>  [15] magrittr_2.0.3              plotly_4.10.4.9000         
#>  [17] sass_0.4.9                  rmarkdown_2.27             
#>  [19] jquerylib_0.1.4             yaml_2.3.8                 
#>  [21] remotes_2.5.0               httpuv_1.6.15              
#>  [23] sctransform_0.4.1           askpass_1.2.0              
#>  [25] spam_2.10-0                 spatstat.sparse_3.0-3      
#>  [27] reticulate_1.37.0           cowplot_1.1.3              
#>  [29] pbapply_1.7-2               RColorBrewer_1.1-3         
#>  [31] abind_1.4-5                 zlibbioc_1.48.0            
#>  [33] Rtsne_0.17                  GenomicRanges_1.54.1       
#>  [35] purrr_1.0.2                 R.utils_2.12.3             
#>  [37] BiocGenerics_0.48.1         RCurl_1.98-1.14            
#>  [39] GenomeInfoDbData_1.2.11     IRanges_2.36.0             
#>  [41] S4Vectors_0.40.2            ggrepel_0.9.5              
#>  [43] irlba_2.3.5.1               listenv_0.9.1              
#>  [45] spatstat.utils_3.1-0        umap_0.2.10.0              
#>  [47] goftest_1.2-3               RSpectra_0.16-1            
#>  [49] spatstat.random_3.2-3       fitdistrplus_1.1-11        
#>  [51] parallelly_1.37.1           leiden_0.4.3.1             
#>  [53] codetools_0.2-20            DelayedArray_0.28.0        
#>  [55] tidyselect_1.2.1            farver_2.1.2               
#>  [57] stats4_4.3.2                spatstat.explore_3.2-7     
#>  [59] jsonlite_1.8.8              BiocNeighbors_1.20.0       
#>  [61] progressr_0.14.0            ggridges_0.5.6             
#>  [63] survival_3.6-4              systemfonts_1.1.0          
#>  [65] dbscan_1.2-0                tools_4.3.2                
#>  [67] ragg_1.3.2                  ica_1.0-3                  
#>  [69] glue_1.7.0                  gridExtra_2.3              
#>  [71] SparseArray_1.2.2           xfun_0.44                  
#>  [73] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.1        
#>  [75] withr_3.0.0                 BiocManager_1.30.23        
#>  [77] fastmap_1.2.0               fansi_1.0.6                
#>  [79] openssl_2.2.0               digest_0.6.35              
#>  [81] rsvd_1.0.5                  R6_2.5.1                   
#>  [83] mime_0.12                   textshaping_0.3.7          
#>  [85] colorspace_2.1-0            scattermore_1.2            
#>  [87] sccore_1.0.5                tensor_1.5                 
#>  [89] spatstat.data_3.0-4         R.methodsS3_1.8.2          
#>  [91] utf8_1.2.4                  tidyr_1.3.1                
#>  [93] generics_0.1.3              data.table_1.15.4          
#>  [95] httr_1.4.7                  htmlwidgets_1.6.4          
#>  [97] S4Arrays_1.2.0              uwot_0.2.2                 
#>  [99] pkgconfig_2.0.3             gtable_0.3.5               
#> [101] lmtest_0.9-40               SingleCellExperiment_1.24.0
#> [103] XVector_0.42.0              htmltools_0.5.8.1          
#> [105] dotCall64_1.1-1             Biobase_2.62.0             
#> [107] png_0.1-8                   SpatialExperiment_1.12.0   
#> [109] knitr_1.46                  reshape2_1.4.4             
#> [111] rjson_0.2.21                nlme_3.1-164               
#> [113] cachem_1.1.0                zoo_1.8-12                 
#> [115] KernSmooth_2.23-24          parallel_4.3.2             
#> [117] miniUI_0.1.1.1              pillar_1.9.0               
#> [119] grid_4.3.2                  vctrs_0.6.5                
#> [121] RANN_2.6.1                  lsa_0.73.3                 
#> [123] promises_1.3.0              xtable_1.8-4               
#> [125] cluster_2.1.6               evaluate_0.23              
#> [127] magick_2.8.3                cli_3.6.2                  
#> [129] compiler_4.3.2              rlang_1.1.3                
#> [131] crayon_1.5.2                future.apply_1.11.2        
#> [133] labeling_0.4.3              mclust_6.1.1               
#> [135] plyr_1.8.9                  stringi_1.8.4              
#> [137] BiocParallel_1.36.0         viridisLite_0.4.2          
#> [139] deldir_2.0-4                munsell_0.5.1              
#> [141] lazyeval_0.2.2              spatstat.geom_3.2-9        
#> [143] Matrix_1.6-5                RcppHNSW_0.6.0             
#> [145] patchwork_1.2.0             future_1.33.2              
#> [147] shiny_1.8.1.1               highr_0.10                 
#> [149] SummarizedExperiment_1.32.0 ROCR_1.0-11                
#> [151] leidenAlg_1.1.3             igraph_2.0.3               
#> [153] bslib_0.7.0
```
</details>
<br>
<br>


