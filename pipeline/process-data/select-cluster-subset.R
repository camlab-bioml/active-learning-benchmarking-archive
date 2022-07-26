library(tidyverse)
library(splitstackshape)
library(SingleCellExperiment)

set.seed(42)

seu <- read_tsv(snakemake@input[['seurat']])
sce <- readRDS(snakemake@input[['sce']])

pred_cell_types <- seu$predicted_cell_type %>% 
  unique() %>% 
  length()

no_to_select <- floor(as.integer(snakemake@wildcards[['cell_num']]) / pred_cell_types)

subset <- stratified(seu, "predicted_cell_type", no_to_select) %>% 
  pull(cell_id)

sce_subset <- sce[, colnames(sce) %in% subset]

if(as.numeric(snakemake@wildcards[['corrupt']]) != 0){
    num_corrupt <- round(as.integer(snakemake@wildcards[['cell_num']]) * as.numeric(snakemake@wildcards[['corrupt']]), 0)
    corrupt_idx <- sample(1:ncol(sce_subset), num_corrupt)
    all_cell_types <- unique(sce$CellType)

    #sce_subset[,corrupt_idx]$CellType <- sample(sce_subset[,corrupt_idx]$CellType)

    for(i in 1:length(all_cell_types)){
        tc <- all_cell_types[i]

        to_corrupt <- which(sce_subset$CellType == tc)
        to_corrupt <- intersect(to_corrupt, corrupt_idx)
        if(length(to_corrupt) > 0){
            subset_cell_types <- all_cell_types[!grepl(tc, all_cell_types)]
            sce_subset$CellType[to_corrupt] <- sample(subset_cell_types, length(to_corrupt), replace = TRUE)
        }
    }
}

gt <- tibble(cell_id = colnames(sce_subset),
             cell_type = sce_subset$CellType,
             method = "Seurat-clustering",
             params = paste0("knn-", snakemake@wildcards[['neighbors']], "-res-", 
                             snakemake@wildcards[['res']], 'cell-num-', snakemake@wildcards[['cell_num']], '-corruption_percentage-', snakemake@wildcards[['corrupt']]))


saveRDS(sce_subset, snakemake@output[['sce']])
write_tsv(gt, snakemake@output[['ground_truth']])
