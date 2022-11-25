library(tidyverse)
library(splitstackshape)
library(SingleCellExperiment)

set.seed(42)

seu <- read_tsv(snakemake@input[['seurat']])
sce <- readRDS(snakemake@input[['sce']])

req_cells <- as.integer(snakemake@wildcards[['cell_num']])

n_cells <- 0

while(n_cells < req_cells){
  pred_cell_types <- seu$predicted_cell_type %>% 
    unique() %>% 
    length()
  
  if(exists("sce_subset")){
    no_to_select <- ceiling((req_cells - ncol(sce_subset)) / pred_cell_types)
  }else{
    no_to_select <- ceiling(req_cells / pred_cell_types)
  }
  
  subset <- stratified(seu, "predicted_cell_type", no_to_select) |> 
    pull(cell_id)
  
  n_cells <- n_cells + length(subset)
  oversampled <- n_cells - req_cells
  if(oversampled > 0){
    subset <- subset[-sample(1:length(subset), oversampled)]
  }
  
  if(exists("sce_subset")){
    sce_subset <- cbind(sce_subset, sce[, colnames(sce) %in% subset])
  }else{
    sce_subset <- sce[, colnames(sce) %in% subset]
  }
  seu <- filter(seu, !(cell_id %in% colnames(sce_subset)))
}

gt <- tibble(cell_id = colnames(sce_subset),
             cell_type = sce_subset$CellType,
             method = "Seurat-clustering",
             params = paste0("knn-", snakemake@wildcards[['neighbors']], "-res-", 
                             snakemake@wildcards[['res']], '-cell-num-', 
                             snakemake@wildcards[['cell_num']]))


saveRDS(sce_subset, snakemake@output[['sce']])
write_tsv(gt, snakemake@output[['ground_truth']])
