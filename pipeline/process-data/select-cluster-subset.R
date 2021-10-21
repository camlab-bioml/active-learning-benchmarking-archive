library(tidyverse)
library(splitstackshape)
library(SingleCellExperiment)

set.seed(42)

seu <- read_tsv(snakemake@input[['seurat']])
sce <- readRDS(snakemake@input[['sce']])

pred_cell_types <- seu$predicted_cell_type %>% 
  unique() %>% 
  length()

no_to_select <- floor(500 / pred_cell_types)

subset <- stratified(seu, "predicted_cell_type", no_to_select) %>% 
  pull(cell_id)

sce_subset <- sce[, colnames(sce) %in% subset]
gt <- tibble(cell_id = colnames(sce_subset),
             cell_type = sce_subset$CellType,
             method = "Seurat-clustering",
             params = paste0("knn-", snakemake@wildcards[['neighbors']], "-res-", snakemake@wildcards[['res']]))


saveRDS(sce_subset, snakemake@output[['sce']])
write_tsv(gt, snakemake@output[['ground_truth']])
