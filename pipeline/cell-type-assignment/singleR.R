suppressPackageStartupMessages({
  library(SingleR)
  library(tidyverse)
  library(caret)
})

train_sce <- readRDS(snakemake@input[['train_data']])
assays(train_sce)$counts <- NULL

train_labels <- read_tsv(snakemake@input[['annotation']]) %>% 
  as.data.frame() %>% 
  column_to_rownames('cell_id')

train_sce$cell_type <- train_labels[colnames(train_sce), 'cell_type']

train_sce <- train_sce[,!is.na(train_sce$cell_type)]

train_sce <- train_sce[,colnames(train_sce) %in% rownames(train_labels)]

annotate_sce <- readRDS(snakemake@input[['test_data']])
assays(annotate_sce)$counts <- NULL

pred <- SingleR(annotate_sce, train_sce, labels = train_sce$cell_type)

result <- tibble(cell_id = rownames(pred),
                 predicted_cell_type = pred$labels,
                 prediction_params = paste0('singleR-labels-iterations_set-', 
                                            snakemake@wildcards[['set']],
                                            "-knn-", snakemake@wildcards[['neighbors']],
                                            "-res-", snakemake@wildcards[['res']], 
                                            "-cell_numbers-", snakemake@wildcards[['cell_num']],
                                            '-randomSelection-', snakemake@wildcards[['rand']], 
                                            '-corrupted-', snakemake@wildcards[['corrupt']]),
                 selection_procedure = paste0(snakemake@wildcards[['selection_procedure']],
                                              '-strategy-', snakemake@wildcards[['strat']]),
                 training_annotator = snakemake@wildcards[['annotator']],
                 modality = 'scRNASeq')

if(is.null(snakemake@wildcards[['cell_selection']])){
  result$cell_selection <- NA
}else{
  result$cell_selection <- snakemake@wildcards[['cell_selection']]
}

write_tsv(result, snakemake@output[['predictions']])