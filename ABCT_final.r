# copy of ABCT_250121.r

# -----------------------------------------------------------
# utils
# -----------------------------------------------------------

#' Filter Duplicated Genes Based on Average Log2 Fold Change
#' 
#' This function filters out duplicated gene entries in a data frame, retaining only the row with 
#' the highest average log2 fold change (avg_log2FC) for each unique gene. The rows are first sorted 
#' in descending order by avg_log2FC, and then duplicates are removed based on the gene column.
#' 
#' @param df A data frame containing gene expression data with a column for `gene` and `avg_log2FC` 
#' representing the gene names and their associated average log2 fold changes, respectively.
#' 
#' @return A data frame with duplicates removed, retaining only the highest avg_log2FC for each gene.
#' 
#' @examples
#' filtered_data <- filter_duplicates(df = gene_expression_data)
#' 
filter_duplicates <- function(df) {
    df <- df[order(-df$avg_log2FC), ]  # Sort by avg_log2FC in descending order
    df <- df[!duplicated(df$gene), ]   # Remove duplicates based on gene column
    return(df)
}



#' Compare Cell Type Annotation Methods Using Accuracy and Weighted F1 Score
#' 
#' This function compares the performance of different cell type annotation methods by calculating 
#' the accuracy and weighted F1 score between the predicted annotations and the actual ground truth labels. 
#' It iterates through a list of predicted annotations and calculates these metrics for each, 
#' returning a summary of the results in a data frame.
#' 
#' @param data A data frame containing both actual ground truth labels and predicted annotations for cell types. 
#' The actual ground truth labels should be in a column specified by the `actual` argument, and the predicted labels
#' should be in columns specified by the `predicted` argument.
#' @param actual The column name in the data containing the actual (ground truth) labels. Default is NULL.
#' @param predicted A list of column names in the data representing the predicted annotations from different methods. Default is NULL.
#' @param digits The number of decimal places to round the accuracy and F1 score results. Default is 2.
#' 
#' @return A data frame summarizing the accuracy and F1 score for each predicted annotation method.
#' 
#' @examples
#' predicted_list <- c("Method1", "Method2", "Method3")
#' result_summary <- compareResult(data = cell_data, actual = "true_labels", predicted = predicted_list, digits = 2)
#' 
compareResult <-function(data,actual=NULL,predicted=NULL,digits = 2) {
    library(mclust)
    library(aricode)
    library(caret)
    library(pROC)
    
    data <- data.frame(lapply(data, as.character), stringsAsFactors = FALSE)

    compare_res <- data.frame(method = character(), acc = numeric(), f1 = numeric())
    
    for (pred in predicted) {
        message(paste("Calculating Accuracy and Metrics for:", pred))
        
        # calculate accuracy, weighted f1
        acc <- calculate_acc(data[[actual]], data[[pred]])
        f1 <- calculate_f1(data[[actual]], data[[pred]])

        compare_res <- rbind(compare_res, data.frame(method = pred, acc = round(acc, digits),  f1 = round(f1, digits)))
    }
    
    print(compare_res)
    return(compare_res)
}


# Calulate accuracy
calculate_acc <- function(actual, predicted) {
  mean(actual == predicted)
}


# Calulate weighted f1 score
calculate_f1 <- function(actual, predicted) {

    all_classes <- union(unique(actual), unique(predicted))

    # Confusion matrix 
    conf_matrix <- table(factor(actual, levels = all_classes),
                         factor(predicted, levels = all_classes))

    # Calculate metrics for each cell type
    tp <- diag(conf_matrix)  # True Positives for each type
    fp <- colSums(conf_matrix) - tp  # False Positives for each type
    fn <- rowSums(conf_matrix) - tp  # False Negatives for each type

    # Avoid division by zero
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)

    # Calculate Weighted Recall
    proportions <- rowSums(conf_matrix) / sum(conf_matrix)
    weighted_recall <- sum(recall * proportions, na.rm = TRUE)

    # Calculate Weighted Precision
    weighted_precision <- sum(precision * proportions, na.rm = TRUE)

    # Calculate Weighted F1 Score
    if (weighted_precision + weighted_recall > 0) {
        weighted_f1 <- 2 * (weighted_precision * weighted_recall) / (weighted_precision + weighted_recall)
    } else {
        weighted_f1 <- 0
    }

    return(weighted_f1)
}



#' Select principal components based on explained variance
#'
#' This function identifies the principal components (PCs) to retain by evaluating 
#' the percentage of variation explained by each PC and selecting the ones that meet 
#' specific variance thresholds. The function considers PCs that explain more than 
#' 90% of cumulative variance and where the variation associated with each PC is less 
#' than a specified threshold.
#'
#' @param data A Seurat object or other object containing dimensional reduction results (e.g., PCA).
#' @param reduc The name of the dimensional reduction technique to evaluate (e.g., "pca"). Default is "pca".
#' @param var The minimum change in variance between consecutive PCs used to identify the significant PCs. Default is 0.05.
#'
#' @return A vector containing two values:
#' \enumerate{
#'   \item The index of the first PC where cumulative variance exceeds 90% and the percentage variance explained by the PC is less than 5%.
#'   \item The index of the last PC where the difference in explained variance between consecutive PCs is greater than the specified threshold (`var`).
#' }
#'
#' @examples
#' pcs <- getPCs(data, reduc = "pca", var = 0.05)
#' print(pcs)
#'
#' @name getPCs
NULL

getPCs <- function(data, reduc = "pca", var = 0.05) {
    # Determine percent of variation associated with each PC
    pct <- data[[reduc]]@stdev / sum(data[[reduc]]@stdev)*100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % 
    # variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    # Last point where change of % of variation is more than 0.05%
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > var), 
                decreasing = T)[1] + 1
    c(co1, co2)
}



#' Create an Area-Assay Based on Metadata and Existing Assay Data
#' 
#' This function creates a new assay called "AreaAssay" by combining the existing assay data (e.g., RNA counts or data) 
#' with metadata information representing a specified "area" from the `meta.data` slot. The `area` information is 
#' replicated according to a specified weight, and the resulting matrix is added as a new assay to the Seurat object. 
#' The function allows the user to specify which layer of data (e.g., "counts" or "data") to incorporate from the existing assay.
#' 
#' @param data A Seurat object containing the assay data and metadata.
#' @param assay The name of the assay from which data will be extracted (default is "RNA").
#' @param layer The type of assay data to use (either "counts" or "data", default is "counts").
#' @param area The column name in the `meta.data` slot that contains the area information to be added as a new assay.
#' @param area_weight The number of times to replicate the area information (default is 1).
#' 
#' @return A Seurat object with the new "AreaAssay" added, incorporating the area information and original assay data.
#' 
#' @examples
#' seurat_obj <- CreateAreaAssay(data = seurat_obj, assay = "RNA", layer = "counts", area = "area_column", area_weight = 1)
#' 
CreateAreaAssay <- function(data, assay = "RNA", layer = "counts", area = NULL, area_weight = 1){
    mat <- GetAssayData(data, assay = assay, layer = layer)
    area_mat <- do.call(rbind, replicate(area_weight, data@meta.data[, area], simplify = FALSE))
	rownames(area_mat) <- paste0("Area.", 1:area_weight)
    mat <- rbind(mat, area_mat)

    if(layer == "counts"){data[["AreaAssay"]] <- CreateAssayObject(counts = mat)}
        else if(layer == "data"){data[["AreaAssay"]] <- CreateAssayObject(data = mat)}

    DefaultAssay(data) <- "AreaAssay"

    return(data)
}



#' Update Metadata in Seurat Object with Additional Information from Another Seurat Object
#' 
#' This function updates the metadata of a Seurat object (`data`) by merging it with metadata from a second Seurat object (`subdata`).
#' It selects specific columns (e.g., `max_score`, `merged_anchor`, `ABCT_clust`, `ABCT_prob`, `ABCT_cut`) from `subdata`'s metadata 
#' and merges them into `data`'s metadata. Missing values (NA) in the `ABCT_clust` and `ABCT_cut` columns are replaced with "Malignant."
#' Optionally, it updates the factor levels of the `ABCT_clust` and `ABCT_cut` columns based on a provided list of cell types (`celltype_list`).
#' 
#' @param data A Seurat object whose metadata will be updated.
#' @param subdata A Seurat object that contains the metadata to merge into `data`.
#' @param celltype_list An optional vector of cell types used to filter and update the factor levels for `ABCT_clust` and `ABCT_cut`.
#' 
#' @return The updated Seurat object (`data`) with the merged and modified metadata.
#' 
#' @examples
#' updated_data <- update_metadata(data = seurat_obj, subdata = sub_seurat_obj, celltype_list = c("T_cells", "B_cells", "Malignant"))
#' 
update_metadata <- function(data, subdata, celltype_list) {

    data$tmp_ID <- rownames(data@meta.data)
    subdata$tmp_ID <- rownames(subdata@meta.data)

    # Select columns from subdata metadata
    tmp <- subdata@meta.data[, c("tmp_ID", "max_score", "merged_anchor", "ABCT_clust", "ABCT_prob", "ABCT_cut")]

    # Merge metadata from data and subdata
    tmp <- left_join(data@meta.data, tmp, by = "tmp_ID")

    # Replace NA values with "Malignant"
    tmp$ABCT_clust <- as.character(tmp$ABCT_clust)
    tmp$ABCT_clust[is.na(tmp$ABCT_clust)] <- "Malignant"

    tmp$ABCT_cut <- as.character(tmp$ABCT_cut)
    tmp$ABCT_cut[is.na(tmp$ABCT_cut)] <- "Malignant"

    if(!is.null(celltype_list)){
        # Update ABCT_clust factor levels and colors
        clust_levelz <- celltype_list[celltype_list %in% tmp$ABCT_clust]
        tmp$ABCT_clust <- factor(tmp$ABCT_clust, levels = clust_levelz)

        # Update ABCT_cut factor levels and colors
        cut_levelz <- celltype_list[celltype_list %in% tmp$ABCT_cut]
        tmp$ABCT_cut <- factor(tmp$ABCT_cut, levels = cut_levelz)
    }

    # Update rownames and assign back to data metadata
    rownames(tmp) <- tmp$tmp_ID; tmp$tmp_ID <- NULL
    data@meta.data <- tmp

    return(data)
}





# -----------------------------------------------------------
# ABCT
# -----------------------------------------------------------
#' Identify and classify anchor cells based on marker gene expression for cell type annotation
#'
#' This function performs supervised classification and identification of anchor cells, 
#' using a gene marker list and expression profiles. It also refines the anchors based 
#' on scoring, and provides options for score smoothing and visualization. The final 
#' output includes malignant anchors, scores, and classification results.
#'
#' @param data A Seurat object containing cell expression data.
#' @param assay The assay to be used from the Seurat object. Default is the active assay.
#' @param ctrl_assay An optional name of the assay representing background counts (e.g., negative controls). Default is NULL.
#' @param marker_list A data frame containing marker genes for cell type classification.
#' @param color_list An optional color list for plotting cell types. Default is NULL.
#' @param cutoff The probability cutoff for classifying cells as "Unknown." Cells with a probability below this value will be labeled as Unknown. Default is NULL.
#' @param method The method used to identify top anchors. Options are "bimodal" or "quantile." Default is "quantile."
#' @param anchor_quantile The quantile threshold used for anchor selection when `method = "quantile"`. Default is 20.
#' @param min_n_anchor The minimum number of anchor cells required for a cell type. Default is NULL.
#' @param use_spatial A logical value indicating whether spatial data should be used. Default is TRUE.
#' @param lambda A parameter for spatial analysis regularization. Default is 0.1.
#' @param k_geom A parameter controlling the neighborhood size for spatial analysis. Default is 15.
#' @param M A parameter influencing the number of features in spatial analysis. Default is 1.
#' @param smooth_k The number of nearest neighbors used for smoothing. Default is 10.
#' @param smooth_reduction The dimensional reduction method used for score smoothing. Default is NULL.
#' @param output_plots Logical value indicating whether to output diagnostic plots. Default is TRUE.
#' @param path The file path where results and plots will be saved. Default is the current working directory.
#'
#' @return A Seurat object with additional metadata, including:
#' \enumerate{
#'   \item max_score: A vector representing the maximum score for each cell based on marker expression.
#'   \item anchor: A factor assigning each cell to a "filter" or anchor category.
#'   \item merged_anchor: The final anchor classification after merging and filtering.
#'   \item ABCT_clust: The cluster assignment for each cell.
#'   \item ABCT_prob: The probability of each cell's cluster assignment.
#'   \item ABCT_cut: The final anchor assignment, including "Unknown."
#' }
#'
#' @examples
#' data <- RunABCT(
#'   data = data, 
#'   assay = "SCT", 
#'   ctrl_assay = "ControlAssay", 
#'   marker_list = markers, 
#'   cutoff = 0.95, 
#'   method = "quantile", 
#'   anchor_quantile = 20, 
#'   min_n_anchor = 5, 
#'   use_spatial = TRUE, 
#'   lambda = 0.1, 
#'   k_geom = 15, 
#'   M = 1, 
#'   smooth_k = 10, 
#'   smooth_reduction = "spatial_pca", 
#'   output_plots = TRUE, 
#'   path = "./results"
#' )
#'
#' @seealso \code{\link{RunBanksy}} , \code{\link{AddModuleScore_UCell}} and \code{\link{SmoothKNN}} for additional functions used in the process.
#'
RunABCT <- function(data, assay = DefaultAssay(data), ctrl_assay = NULL, marker_list = NULL, color_list = NULL, cutoff = NULL, method = "quantile", anchor_quantile = 20, min_n_anchor = NULL, use_spatial = TRUE, lambda = 0.1, k_geom=15, M = 1, w_neg=NULL, dimx=NULL, dimy=NULL, smooth_k = 10, smooth_reduction = NULL, output_plots = TRUE, path = getwd()){

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
    library(Banksy)

    DefaultAssay(data) <- assay

    marker_list$cluster <- factor(marker_list$cluster, levels = unique(marker_list$cluster))

    celltype_list <- names(table(marker_list$cluster))

    #### ERROR ---------------------------
    if(is.null(marker_list)){stop("Error: No marker list")}
    if(length(data@reductions) == 0 & (!is.null(smooth_reduction) | output_plots)){stop("Error: No reduction")}
    if(length(color_list) > 0 & (length(color_list) != length(celltype_list))){stop("Error: Different lengths of celltype_list & color_list")}
    if(output_plots & !is.null(color_list) & any(!grepl("^#", names(color_list)))){stop("Error: Invalid color code")}


    #### Setting plot options ---------------------------
    if(output_plots){
        ifelse(use_spatial, banksy_subTitle <- paste("BANKSY k_geom = ", k_geom, ", lambda = ", lambda, ", M = ", M), banksy_subTitle <- "BANKSY X")
        ifelse(method == "bimodal", method_subTitle <- "method = bimodal", method_subTitle <- paste0("method = quantile ", anchor_quantile, "%"))
        
        ifelse(!is.null(color_list), names(celltype_list) <- color_list, names(celltype_list) <- hue_pal()(length(celltype_list)))
    }

    ### Creating SpatialAssay --------------------------
    if(use_spatial & !("BANKSY" %in% names(data@assays))){
        message("Creating SpatialAssay")
        if(is.data.frame(M)){
            max.M <- max(M$m, na.rm = T)
            data <- RunBanksy(data, assay = assay, slot = 'data', lambda = lambda, dimx = dimx, dimy = dimy, features = 'all', k_geom = k_geom, M = max.M)
        } else if(is.numeric(M)){
            data <- RunBanksy(data, assay = assay, slot = 'data', lambda = lambda, dimx = dimx, dimy = dimy, features = 'all', k_geom = k_geom, M = M)
        }

        if(smooth_reduction %in% c("spatial_pca", "spatial_umap")){
            data <- RunPCA(data, assay = 'BANKSY', features = rownames(data), npcs = 50, reduction.name = 'spatial_pca', verbose = FALSE)

            if(smooth_reduction == "spatial_umap"){
                npcs <- min(getPCs(data, reduc = "spatial_pca"))
                data <- RunUMAP(data, reduction = 'spatial_pca', reduction.name = 'spatial_umap', dims = 1:npcs, verbose = FALSE, umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
            }
        }

    }

    if(use_spatial){
        DefaultAssay(data) <- "BANKSY"
    }

    ### Calculating UCell module score --------------------------
    message("Calculating UCell module score")
    if(use_spatial){
        if(is.numeric(M)){
            res <- marker_list
            tmp <- marker_list
            tmp$gene <- paste0(tmp$gene, ".m", 0)
            res <- rbind(res, tmp)

            marker_list <- res; rm(res, tmp)
        } else if(is.data.frame(M)){
            res <- marker_list

            for(cellType in unique(M$celltype)){
                tmp.m <- M[M$celltype == cellType, "m"]
                tmp <- marker_list %>% filter(cluster == cellType)
                tmp <- data.frame(
                    "cluster" = cellType,
                    "gene" = paste0(rep(tmp$gene, each = tmp.m + 1), ".m", 0:tmp.m)
                )
                res <- rbind(res, tmp)
            }

            marker_list <- res; rm(res, tmp, tmp.m)
        }
    }
    
    #### Calculate AddModuleScore_UCell ---------------------------
    marker_list <- split(marker_list$gene, marker_list$cluster)

    ## message
    for (cluster in names(marker_list)) {
        genes <- paste(marker_list[[cluster]], collapse = ", ")
        message(cluster, ": ", "\n", genes, "\n")
    }

    data <- AddModuleScore_UCell(data, features = marker_list, slot="data", name = "_score", w_neg=w_neg)

    DefaultAssay(data) <- assay

    #### Max score ---------------------------
    if(!is.null(smooth_reduction)){
        data <- SmoothKNN(data, signature.names = paste0(celltype_list, "_score"), reduction = smooth_reduction, k = smooth_k)
        tmp <- colnames(data@meta.data)[grepl("_kNN", colnames(data@meta.data))]
    } else{
        tmp <- colnames(data@meta.data)[grepl("_score", colnames(data@meta.data))]
    }

    find_max_column <- function(row) {
        max_column <- which.max(row)
        return(names(row)[max_column])
    }

    scr_result <- data@meta.data[, tmp]
    max_score <- apply(scr_result, 1, find_max_column)
    max_score<- gsub("_score","",max_score) 
    max_score<- gsub("_kNN","",max_score)

    data$max_score <- max_score

    celltype_list <- celltype_list[celltype_list %in% names(table(data$max_score))]
    data$max_score <- factor(data$max_score, levels = celltype_list)
    message("** CellType max score result **")
    print(table(data$max_score, useNA = "always"))

    if(output_plots){
    	## DimPlot of max score cell type
    	DimPlot(data, raster=FALSE, label = TRUE, group.by= "max_score", label.size=5, pt.size = 0.001, reduction = "umap", cols = names(celltype_list)) +
            ggtitle(paste0("Max score (",banksy_subTitle, ")"))
    	ggsave(paste0(path, "/1.Dimplot_celltype_max_score.png"), width = 12, height = 8)        
    }

    #### Find top anchors ---------------------------
    if(method == "bimodal"){
        for(celltype in celltype_list){ 
            tmp <- data@meta.data %>% filter(max_score == celltype)
            
            if(!is.null(smooth_reduction)){
                module_score <- tmp[[paste0(celltype, "_score", "_kNN")]]
            } else{
                module_score <- tmp[[paste0(celltype, "_score")]]
            }

            km <- kmeans(module_score, centers = 2)
            d <- density(module_score)
            v <- optimize(approxfun(d$x,d$y),interval = km$centers[,1])$minimum

            # Draw malignant score density plot (나중에 뺄 수 있는거...)
            ggplot(data = data.frame(setNames(list(module_score), paste0(celltype, "_score"))), aes(x = get(paste0(celltype, "_score")))) +
            geom_density() +
                geom_vline(xintercept = km$centers[,1], color = "blue") +
                labs(x = paste0(celltype, "_score"), y = "density") +
                ggtitle(paste0("merged_anchor (", method_subTitle, ")"))
            ggsave(paste0(path,"/2.",celltype,"_score_densityPlot.png"), width = 10)

            celltype_anchor <- rep("filter", dim(tmp)[1])
            names(celltype_anchor) <- rownames(tmp)

            celltype_anchor[module_score > max(km$centers[,1])] <- celltype

            anchor <- rep("filter", dim(data)[2])
            names(anchor) <- colnames(data)
        
            anchor[rownames(tmp)] <- celltype_anchor[rownames(tmp)]

            message(paste0("* ", celltype, " anchor"))
            print(table(anchor))

            anchor <- factor(anchor, levels = c(celltype, "filter"))
            data@meta.data[[paste0(celltype, "_anchor")]] <- anchor

            if(output_plots){ ## DimPlot of anchors
                DimPlot(data,raster=FALSE, label = TRUE, group.by= paste0(celltype, "_anchor"),label.size=5, pt.size = 0.001, reduction = "umap", cols = c(names(celltype_list)[celltype_list == celltype], "gray")) +
                ggtitle(paste0(celltype, "_anchor (", banksy_subTitle, ", ", method_subTitle, ")"))
                ggsave(paste0(path, "/2.Dimplot_", celltype,"_anchor.png"), width = 12, height = 8)
            }
	    } 

    } else if (method == "quantile"){
        for(celltype in celltype_list){ 
            tmp <- data@meta.data %>% filter(max_score == celltype)
            
            if(!is.null(smooth_reduction)){
                module_score <- tmp[[paste0(celltype, "_score", "_kNN")]]
            } else{
                module_score <- tmp[[paste0(celltype, "_score")]]
            }

        	celltype_anchor <- rep("filter", dim(tmp)[1])
        	names(celltype_anchor) <- rownames(tmp)

        	celltype_anchor[module_score >= quantile(module_score, (1-anchor_quantile/100))] <- celltype

        	anchor <- rep("filter", dim(data)[2])
        	names(anchor) <- colnames(data)
        
        	anchor[rownames(tmp)] <- celltype_anchor

        	message(paste0("* ", celltype, " anchor"))
            print(table(anchor))

        	anchor <- factor(anchor, levels = c(celltype, "filter"))
        	data@meta.data[[paste0(celltype, "_anchor")]] <- anchor

        	if(output_plots){## DimPlot of anchors
                ifelse(length(unique(data@meta.data[[paste0(celltype, "_anchor")]])) == 1, colz <- "gray", colz <- c(names(celltype_list)[celltype_list == celltype], "gray"))

                DimPlot(data,raster=FALSE, label = TRUE, group.by= paste0(celltype, "_anchor"),label.size=5, pt.size = 0.001, reduction = "umap", cols = colz) +
                ggtitle(paste0(celltype, "_anchor (", banksy_subTitle, ", ", method_subTitle, ")"))
                ggsave(paste0(path, "/2.Dimplot_", celltype,"_anchor.png"), width = 12, height = 8)
            }
        }
    }
    
    tmp <- colnames(data@meta.data)[grep(paste(celltype_list, collapse = "|"), colnames(data@meta.data))]
    anchor_result <- data@meta.data[, tmp]

    #### Merging anchor result ---------------------------
    anchor_df <- anchor_result[, grep("_anchor", colnames(anchor_result))]

    count_non_filter <- function(row){sum(row != "filter")}

    tmp <- apply(anchor_df, 1, count_non_filter)
    tmp <- (tmp > 1)
    anchor_df[tmp, ] <- "filter"

    merged_anchor <- apply(anchor_df, 1, function(row){
        if(all(row == "filter")) {
            return("filter")
        } else {
            return(row[row != "filter"][1])
        }
    })

    if(!is.null(min_n_anchor)){
        message("Anchor before filtering")
        print(table(merged_anchor, useNA = "always"))

        tmp <- names(table(merged_anchor)[table(merged_anchor) < min_n_anchor])
        celltype_list <- celltype_list[!(celltype_list %in% tmp)]

        for(filt in tmp){
            merged_anchor[merged_anchor == filt] <- "filter"
        }

        message(paste0("Final anchor (count >= ", min_n_anchor,")"))
        print(table(merged_anchor, useNA = "always"))

        message("Filtered cell type")
        message(paste(tmp, collapse = ", "))
    } else {
        message("Final anchor")
        print(table(merged_anchor, useNA = "always"))
    }

    anchor_result <- cbind(anchor_result, merged_anchor)
    write.csv(anchor_result, paste0(path, "/AnchoR_score_result_df.csv"))

    celltype_list <- celltype_list[celltype_list %in% unique(merged_anchor)]
    data$merged_anchor <- factor(merged_anchor, levels = c(celltype_list, "filter"))

    if(output_plots){
        ## DimPlot of anchors
        DimPlot(data, raster=FALSE, label = TRUE, group.by= "merged_anchor", label.size=5, pt.size = 0.001, reduction = "umap", cols = c(names(celltype_list), "gray"), repel = TRUE) +
            ggtitle(paste0("Merged Anchors (", banksy_subTitle, ", ", method_subTitle, ")"))
        ggsave(paste0(path,"/3.Dimplot_merged_anchors.png"), width = 12, height = 8) 
    }

    #### Run insitutype ---------------------------
    message("Running Insitutype")

    if(!is.null(ctrl_assay)){
        neg <- GetAssayData(data, assay = ctrl_assay, layer = "counts") %>% as.matrix
        data$negmean <- Matrix::colMeans(neg)
    }
    else{
        data$negmean <- 0
    }

    ## Create anchor cell profile
    Idents(data) <- as.character(data$merged_anchor)
    cells <- WhichCells(data, ident = c('filter'), invert = TRUE)

    query <- GetAssayData(data, assay = assay, layer = "counts") %>% as.matrix
    counts <- query[, cells]

    clust <- as.character(data$merged_anchor)
    names(clust) <- Cells(data)
    clust <- clust[names(clust) %in% cells]

    neg <- data$negmean
    names(neg) <- Cells(data)
    neg <- neg[names(neg) %in% cells]

    profiles <- getRNAprofiles(x = t(counts), clust = clust, neg = neg)

    anchor_result = insitutypeML(x = t(query),
                        neg = data@meta.data[, "negmean"],
                        assay_type = "RNA", 
                        reference_profiles = profiles,
                        cohort = NULL)
        
    anchor_df <- cbind(anchor_result$clust, anchor_result$prob, anchor_result$clust)
    colnames(anchor_df) <- c("ABCT_clust", "ABCT_prob", "ABCT_cut")
    anchor_df <- as.data.frame(anchor_df)
    anchor_df$ABCT_prob <- as.numeric(anchor_df$ABCT_prob)

    if(is.null(cutoff)){
        res_df <- anchor_df[,c("ABCT_prob"), drop = FALSE]
        res_df$ABCT_prob2 <- round(res_df$ABCT_prob,2)

        # 밀도 계산 (0.5에서 1 사이 구간)
        density_values <- density(res_df$ABCT_prob2, from = 0.5, to = 1, n = 500)  # n은 구간을 0.01 단위로 생성

        # 각 지점에서 밀도의 차이 계산
        density_diff <- diff(density_values$y)  # 인접한 지점들 간의 밀도 변화량 계산

        # 가장 큰 변화가 발생한 지점 찾기
        max_diff_index <- which.max(abs(density_diff))
        max_diff_x <- density_values$x[max_diff_index]
        max_diff_y <- density_values$y[max_diff_index]

        # 결과 출력
        cutoff <- round(max_diff_x,2)
        message(paste0("Probability cutoff: ", cutoff))
    } 
    anchor_df$ABCT_prob2 <- res_df$ABCT_prob2
    anchor_df$ABCT_cut[anchor_df$ABCT_prob2 < cutoff] <- "Unknown"
    anchor_df$ABCT_prob2 <- NULL

    message("Anchor result (clust)")
    print(table(anchor_df$ABCT_clust))

    message("Anchor result (cut)")
    print(table(anchor_df$ABCT_cut))

    message("Saving annotation result as dataframe (csv)")
    write.csv(anchor_df, paste0(path,"AnchoR_insitu_result.csv"))

    message("Saving annotation result as rds object")
    saveRDS(anchor_result, paste0(path,"AnchoR_insitu_result.rds"))

    celltype_list <- celltype_list[celltype_list %in% unique(anchor_df$ABCT_clust)]

    data@meta.data <- cbind(data@meta.data, anchor_df)
    data$ABCT_clust <- factor(data$ABCT_clust, levels = celltype_list)
    data$ABCT_cut <- factor(data$ABCT_cut, levels = c(celltype_list, "Unknown"))

    if(output_plots){
        ## DimPlot
        DimPlot(data, raster=FALSE, label = TRUE, group.by= "ABCT_clust", label.size=5, pt.size = 0.001, reduction = "umap", cols = names(celltype_list), repel = TRUE) +
            ggtitle(paste0("AnchoR final result (", banksy_subTitle, ", ", method_subTitle, ")"))
        ggsave(paste0(path,"/4.Dimplot_ABCT_clust.png"), width = 12, height = 8) 

        DimPlot(data, raster=FALSE, label = TRUE, group.by= "ABCT_cut", label.size=5, pt.size = 0.001, reduction = "umap", cols = c(names(celltype_list), "gray"), repel = TRUE) +
            ggtitle(paste0("AnchoR final result (", banksy_subTitle, ", ", method_subTitle, ")"))
        ggsave(paste0(path,"/4.Dimplot_ABCT_cut.png"), width = 12, height = 8)

        if("spatial_umap" %in% names(data@reductions)){
            ## DimPlot
            DimPlot(data, raster=FALSE, label = TRUE, group.by= "ABCT_clust", label.size=5, pt.size = 0.001, reduction = "spatial_umap", cols = names(celltype_list), repel = TRUE) +
                ggtitle(paste0("AnchoR final result (", banksy_subTitle, ", ", method_subTitle, ")"))
            ggsave(paste0(path,"/4.Dimplot_ABCT_clust_spatial_umap.png"), width = 12, height = 8) 

            DimPlot(data, raster=FALSE, label = TRUE, group.by= "ABCT_cut", label.size=5, pt.size = 0.001, reduction = "spatial_umap", cols = c(names(celltype_list), "gray"), repel = TRUE) +
                ggtitle(paste0("AnchoR final result (", banksy_subTitle, ", ", method_subTitle, ")"))
            ggsave(paste0(path,"/4.Dimplot_ABCT_cut_spatial_umap.png"), width = 12, height = 8)

        }

        # ImgdimPlot
        coord <- GetTissueCoordinates(data)
        ratio <- c(range(coord$x)[2] - range(coord$x)[1], range(coord$y)[2] - range(coord$y)[1])
        max_value <- 14; 
        scale_factor <- max_value / max(ratio); scaled_ratio <- ratio * scale_factor

        ImageDimPlot(data, fov = "fov", group.by="ABCT_clust",
                size = 0.85, dark.background = FALSE, cols = names(celltype_list)) + coord_flip()+
                theme(legend.text=element_text(size=14), legend.title=element_text(size=20)) +
            ggtitle(paste0("AnchoR final result (", banksy_subTitle, ", ", method_subTitle, ")"))
        ggsave(paste0(path, "/5.ImgDim_ABCT_clust.png"), width = (scaled_ratio[1]+2), height = scaled_ratio[2])

        ImageDimPlot(data, fov = "fov", group.by="ABCT_cut",
                size = 0.85, dark.background = FALSE, cols = c(names(celltype_list), "gray")) + coord_flip()+
                theme(legend.text=element_text(size=14), legend.title=element_text(size=20)) +
            ggtitle(paste0("AnchoR final result (", banksy_subTitle, ", ", method_subTitle, ")"))
        ggsave(paste0(path, "/5.ImgDim_ABCT_cut.png"), width = (scaled_ratio[1]+2), height = scaled_ratio[2])

    }


    ### Remove unnecessary metadata columns
    tmp <- c(paste0(names(marker_list), "_anchor"), paste0(names(marker_list), "_score"), paste0(names(marker_list), "_score_kNN"))
    data@meta.data[,tmp] <- NULL

    return(data)

}



#' Identify and classify malignant cells based on gene expression profiles and spatial data
#'
#' This function performs the identification and classification of malignant cells in a Seurat object,
#' using a combination of gene marker lists, spatial information, and score-based methods. The function 
#' allows for the creation of new assays based on area-weighted features, scoring of cells for malignant 
#' potential, and the smoothing of scores for better visualization and clustering. It also provides options 
#' for density-based peak detection to classify malignant and non-malignant cells. 
#' Additionally, the function uses Quadratic Discriminant Analysis (QDA) on PCA embeddings to classify 
#' the cells as malignant or non-malignant. The final output includes the classification of malignant and 
#' non-malignant cells based on the QDA model, along with their scores and clustering results.
#'
#' @param data A Seurat object containing the cell expression data to be analyzed.
#' @param assay The assay to be used for analysis. Default is the active assay.
#' @param ctrl_assay An optional assay name to be used as a control for background counts. Default is NULL.
#' @param marker_list A data frame containing the gene markers for malignant cell classification.
#' @param area An optional value representing the area for defining area-weighted features. Default is NULL.
#' @param area_weight A weighting factor for the area, affecting the influence of area features on the analysis. Default is 1.
#' @param use_spatial A logical value indicating whether to use spatial data. Default is TRUE.
#' @param lambda A regularization parameter for spatial analysis. Default is 0.1.
#' @param k_geom A parameter that defines the neighborhood size for spatial analysis. Default is 15.
#' @param M A parameter influencing the number of features in spatial analysis. Default is 1.
#' @param dimx The spatial x dimension for spatial analysis. Default is NULL.
#' @param dimy The spatial y dimension for spatial analysis. Default is NULL.
#' @param smooth_reduction The dimensionality reduction method used for score smoothing. Default is NULL.
#' @param smooth_k The number of nearest neighbors used for smoothing. Default is 20.
#' @param output_plots A logical value indicating whether to output diagnostic plots. Default is TRUE.
#' @param path The file path where results and plots will be saved. Default is the current working directory.
#' @param w_neg A weight factor for negative scores during score calculation. Default is 1.
#' @param anchor_sample An optional value indicating the number of sampled cells for the anchor classification. Default is NULL.
#'
#' @return A Seurat object with additional metadata, including:
#' \enumerate{
#'   \item malignant_score: A score representing the likelihood of each cell being malignant based on marker expression.
#'   \item malignant_anchor: A factor assigning each cell to "Malignant," "Non_malignant," or "filter" based on its score.
#'   \item malignant_result: The final classification of each cell as "Malignant" or "Non_malignant" based on QDA.
#' }
#'
#' @examples
#' # Example of running FindMalignantCells with various parameters
#' malignant_result <- FindMalignantCells(
#'   data = seurat_data, 
#'   marker_list = markers, 
#'   area = 50, 
#'   area_weight = 2, 
#'   use_spatial = TRUE, 
#'   lambda = 0.1, 
#'   k_geom = 15, 
#'   M = 1, 
#'   smooth_reduction = "spatial_pca", 
#'   smooth_k = 20, 
#'   output_plots = TRUE, 
#'   path = "./results"
#' )
#' 
#' # View the malignant classification results
#' print(table(malignant_result$malignant_result))
#' 
#' @seealso \code{\link{RunBanksy}} and \code{\link{AddModuleScore_UCell}} for additional functions used in the process.
NULL
FindMalignantCells <- function(data, assay = DefaultAssay(data), ctrl_assay=NULL, marker_list = NULL, area = NULL, area_weight = 1, use_spatial = TRUE, lambda = 0.1, k_geom=15, M = 1, dimx=NULL, dimy=NULL, smooth_reduction = NULL, smooth_k = 20, output_plots = TRUE, path = getwd(), w_neg = 1, anchor_sample = NULL){

    #### ERROR ---------------------------
    if(is.null(marker_list)){stop("Error: No marker list")}
    if(length(data@reductions) == 0 & (!is.null(smooth_reduction) | output_plots)){stop("Error: No reduction")}

	library(Seurat)
	library(SeuratWrappers)
	library(harmony)
	library(dplyr)
	library(scales)
	library(EnvStats)
	library(stringr)
	library(matrixStats)
	library(InSituType)
	library(UCell)
    library(Banksy)

    DefaultAssay(data) <- assay
    tmp_slot <- "counts"

    if(output_plots){
        ifelse(!is.null(area), area_subTitle <- paste("Area weight = ", area_weight), area_subTitle <- "Area X")
        ifelse(use_spatial, banksy_subTitle <- paste("BANKSY k_geom = ", k_geom, ", lambda = ", lambda, ", M = ", M), banksy_subTitle <- "BANKSY X")
    }

    #### Creating New Assay ---------------------------
    if(use_spatial & !("BANKSY" %in% names(data@assays))){
        message("Creating SpatialAssay")
        data <- RunBanksy(data, assay = DefaultAssay(data), slot = 'data', lambda = lambda, dimx = dimx, dimy = dimy, features = 'all', k_geom = k_geom, M = M)

        if(!is.null(smooth_reduction)){
            if(smooth_reduction %in% c("spatial_pca", "spatial_umap")){
                data <- RunPCA(data, assay = 'BANKSY', features = rownames(data), npcs = 50, reduction.name = 'spatial_pca', verbose = FALSE)

                if(smooth_reduction == "spatial_umap"){
                    npcs <- min(getPCs(data, reduc = "spatial_pca"))
                    data <- RunUMAP(data, reduction = 'spatial_pca', reduction.name = 'spatial_umap', dims = 1:npcs, verbose = FALSE, umap.method = 'uwot', n.neighbors = 30, n.epochs=-1, min.dist = 0.1)
                }
            }
        }

    }

    if(use_spatial){
        DefaultAssay(data) <- "BANKSY"; tmp_slot <- "data"
    }


    if(!is.null(area)){
        message("Creating AreaAssay")
        data <- CreateAreaAssay(data, assay = DefaultAssay(data), layer = tmp_slot, area = area, area_weight = area_weight)
    }

    #### Defining marker list again ---------------------------
    if(use_spatial){
        res <- marker_list

        tmp <- marker_list
        tmp$gene <- paste0(tmp$gene, ".m", 0)
        res <- rbind(res, tmp)

        marker_list <- res; rm(res, tmp)
    }

    if(!is.null(area)){
        marker_list <- rbind(marker_list, data.frame(cluster = "Malignant", gene = paste0("Area.", 1:area_weight)))
    }

    #### Calculate AddModuleScore_UCell ---------------------------
    marker_list <- split(marker_list$gene, marker_list$cluster)

    ## message
    for (cluster in names(marker_list)) {
        genes <- paste(marker_list[[cluster]], collapse = ", ")
        message(cluster, ": ", "\n", genes, "\n")
    }

    data <- AddModuleScore_UCell(data, features = marker_list, slot = tmp_slot, name = "_score", w_neg = w_neg)
    tmp_names <- paste0(cluster, "_score")
    colnames(data@meta.data)[colnames(data@meta.data) == tmp_names] <- "malignant_score"

    DefaultAssay(data) <- assay

    if(sum(grepl("AreaAssay", names(data@assays)))){
        data[["AreaAssay"]] <- NULL
        data$nCount_AreaAssay <- NULL
        data$nFeature_AreaAssay <- NULL
    }

    #### Score smoothing ---------------------------
    if(!is.null(smooth_reduction)){
        message("Smoothing module score (UCell)")
        data <- SmoothKNN(data, signature.names = "malignant_score", reduction = smooth_reduction, k = smooth_k) # **
        data$malignant_score <- data$malignant_score_kNN
        data$malignant_score_kNN <- NULL
    }

    # ++++++++++++++++++++++++++
    malignant_density <- density(data$malignant_score)
    peaks <- findpeaks(malignant_density$y, threshold = 0.1)

    # peak가 3개 초과하면 density 값 상위 3개만 갖고 오기
    if(nrow(peaks) > 3){
        peaks <- peaks[order(-peaks[, 1]), ]
        peaks <- peaks[1:3,]
    }

    n_peak <- nrow(peaks)
    n_cluster <- ifelse(n_peak <= 2, 2, 3)
    peak_positions <- sort(malignant_density$x[peaks[, 2]])

    ## anchor select - peak  ---------------------------
    if(n_peak == 1){
        #### GMM ---------------------------
        gmm_model <- Mclust(data$malignant_score, G = n_cluster)
        median_values <- tapply(data$malignant_score, gmm_model$classification, median) %>% sort

        ## Find malignant anchors
        malignant_anchor <- rep("filter", dim(data)[2])
        names(malignant_anchor) <- colnames(data)
        malignant_anchor[data$malignant_score > median_values[2]] <- "Malignant"
        malignant_anchor[data$malignant_score < median_values[1]] <- "Non_malignant"
        data$malignant_anchor <- factor(malignant_anchor, levels = c("Malignant", "Non_malignant", "filter"))    
    } else{
        ## Find malignant anchors
        malignant_anchor <- rep("filter", dim(data)[2])
        names(malignant_anchor) <- colnames(data)
        malignant_anchor[data$malignant_score > peak_positions[n_peak]] <- "Malignant"
        malignant_anchor[data$malignant_score < peak_positions[n_peak-1]] <- "Non_malignant"
        data$malignant_anchor <- factor(malignant_anchor, levels = c("Malignant", "Non_malignant", "filter"))        
    }

    # #### train & test data (PCA embedding) ----------------
    if (use_spatial) {
        mclust_dimreduc <- "spatial_pca" 
        if(is.null(smooth_reduction)){
            mclust_dimreduc <- "pca"
        } ## 지우
    } else {
        mclust_dimreduc <- "pca"
    }

    # if (smooth_reduction=="pca"){  #지원
    #     mclust_dimreduc <- "pca"
    # }

    Idents(data) <- as.character(data$malignant_anchor)

    npcs <- min(getPCs(data, reduc = mclust_dimreduc))
    message(paste0("npcs: ", npcs))
    # query <- data@reductions[[mclust_dimreduc]]@cell.embeddings %>% as.data.frame
    query <- data@reductions[[mclust_dimreduc]]@cell.embeddings[, 1:npcs] %>% as.data.frame

    cells <- WhichCells(data, ident = c('filter'), invert = TRUE)
    clust <- as.character(data$malignant_anchor)
    names(clust) <- Cells(data)
    clust <- clust[cells]

    if(!is.null(anchor_sample)){
        clust <- c(sample(clust[clust == "Malignant"], min(anchor_sample, sum(clust == "Malignant"))), sample(clust[clust == "Non_malignant"], min(anchor_sample, sum(clust == "Non_malignant"))))
        clust_tmp <- data.frame(tmp_key = names(clust), malignant_anchor_sampled = clust)
        tmp_meta <- data.frame(tmp_key = rownames(data@meta.data))
        tmp_meta <- left_join(tmp_meta, clust_tmp, by = "tmp_key")
        data$malignant_anchor_sampled <- tmp_meta$malignant_anchor_sampled
        data$malignant_anchor_sampled[is.na(data$malignant_anchor_sampled)] <- "filter"
        data$malignant_anchor_sampled <- factor(data$malignant_anchor_sampled, levels = c("Malignant", "Non_malignant", "filter"))
        rm(clust_tmp, tmp_meta)
    }

    # train_data
    train_data <- query[names(clust), ]
    train_data$clust <- factor(clust, levels = c("Malignant", "Non_malignant"))

    #### * QDA (Quadratic Discriminant Analysis) -------------------------
    library(MASS)

    qda_model <- qda(clust ~ ., data = train_data) # model training
    data$malignant_result <- predict(qda_model, query)$class

    data$malignant_result <- factor(data$malignant_result, levels = c("Malignant", "Non_malignant"))

    if(output_plots){
        # Density Plot - QDA
        density_df <- data.frame(x = malignant_density$x, y = malignant_density$y)

        cluster_density <- data.frame(cluster = factor(data$malignant_result),
                                        score = data$malignant_score) %>%
                        group_by(cluster) %>%
                        do(data.frame(density_x = density(.$score)$x, 
                                        density_y = density(.$score)$y))

        peak_positions_df <- data.frame(x = sort(malignant_density$x[peaks[, 2]]))

        ggplot() +
        geom_line(data = density_df, aes(x = x, y = y), color = "black", size = 1) +
        geom_area(data = cluster_density, aes(x = density_x, y = density_y, fill = cluster), alpha = 0.4, color = NA) +
        geom_line(data = cluster_density, aes(x = density_x, y = density_y, color = cluster), size = 1) +
        geom_vline(data = peak_positions_df, aes(xintercept = x), color = "blue", linetype = "solid", size = 1) +
        labs(title = paste0("Malignant Score Density Plot"),
            x = "Malignant Score",
            y = "Density")
        ggsave(paste0(path,"/1.Malignant_score_densityPlot.png"), width = 10)

        # DimPlot
        if(!is.null(anchor_sample)){
            DimPlot(data,raster=FALSE,label = TRUE, group.by= "malignant_anchor_sampled",label.size=5, pt.size = 0.001, reduction = "umap", cols = c(hue_pal()(2), "gray88"))+
                ggtitle(paste0("Malignant Anchor (Sampled)"))
            ggsave(paste0(path, "/2.Dimplot_malignant_anchor_sampled.png"), width = 12, height = 8)
        }

        DimPlot(data,raster=FALSE,label = TRUE, group.by= "malignant_anchor",label.size=5, pt.size = 0.001, reduction = "umap", cols = c(hue_pal()(2), "gray88"))+
            ggtitle(paste0("Malignant Anchor"))
        ggsave(paste0(path, "/2.Dimplot_malignant_anchor.png"), width = 12, height = 8)

        DimPlot(data,raster=FALSE,label = TRUE, group.by= "malignant_result",label.size=5, pt.size = 0.001, reduction = "umap", cols = hue_pal()(2))+
            ggtitle(paste0("Malignant Result"))
        ggsave(paste0(path, "/2.Dimplot_malignant_result.png"), width = 12, height = 8)

        if("spatial_umap" %in% names(data@reductions)){
            DimPlot(data,raster=FALSE,label = TRUE, group.by= "malignant_anchor",label.size=5, pt.size = 0.001, reduction = "spatial_umap", cols = c(hue_pal()(2), "gray88"))+
                ggtitle(paste0("Malignant Anchor"))
            ggsave(paste0(path, "/2.Dimplot_malignant_anchor_spatial_umap.png"), width = 12, height = 8)

            DimPlot(data,raster=FALSE,label = TRUE, group.by= "malignant_result",label.size=5, pt.size = 0.001, reduction = "spatial_umap", cols = hue_pal()(2))+
                ggtitle(paste0("Malignant Result"))
            ggsave(paste0(path, "/2.Dimplot_malignant_result_spatial_umap.png"), width = 12, height = 8)            
        }

        # ImgdimPlot
        coord <- GetTissueCoordinates(data)
        ratio <- c(range(coord$x)[2] - range(coord$x)[1], range(coord$y)[2] - range(coord$y)[1])
        max_value <- 14; 
        scale_factor <- max_value / max(ratio); scaled_ratio <- ratio * scale_factor

        ImageDimPlot(data, fov = "fov", group.by="malignant_result",
                size = 0.85, dark.background = FALSE) + coord_flip()+
                theme(legend.text=element_text(size=14), legend.title=element_text(size=20)) +
        ggtitle(paste0("Malignant Result"))
        ggsave(paste0(path, "/3.ImgDim_malignant_result.png"), width = (scaled_ratio[1]+2), height = scaled_ratio[2])

    }

    return(data)
}









