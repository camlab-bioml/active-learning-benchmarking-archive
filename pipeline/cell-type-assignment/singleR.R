suppressPackageStartupMessages({
  library(SingleR)
  library(tidyverse)
})

train_sce <- readRDS(snakemake@input[['train_data']])
assays(train_sce)$counts <- NULL

train_labels <- read_tsv(snakemake@input[['annotation']]) %>% 
  as.data.frame() %>% 
  column_to_rownames('cell_id')

train_sce$cell_type <- train_labels[colnames(train_sce), 'cell_type']

train_sce <- train_sce[,!is.na(train_sce$cell_type)]

annotate_sce <- readRDS(snakemake@input[['test_data']])
assays(annotate_sce)$counts <- NULL

pred <- SingleR(annotate_sce, train_sce, labels = train_sce$cell_type)

result <- tibble(cell_id = rownames(pred),
                 predicted_cell_type = pred$labels,
                 prediction_params = paste0('singleR-labels-iterations_set-', 
                                            snakemake@wildcards[['set']],
                                            "-knn-", snakemake@wildcards[['neighbors']],
                                            "-res-", snakemake@wildcards[['res']]),
                 selection_procedure = snakemake@wildcards[['selection_procedure']],
                 training_annotator = snakemake@wildcards[['annotator']],
                 modality = 'scRNASeq')

write_tsv(result, snakemake@output[['predictions']])