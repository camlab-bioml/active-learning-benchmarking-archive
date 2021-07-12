suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(caret)
})

sce <- readRDS(snakemake@input[['rds']])

set.seed(42)
train <- createDataPartition(sce$CellType, p = 0.5)$Resample1

train_sce <- sce[,train]
test_sce <- sce[,-train]

saveRDS(train_sce, snakemake@output[['train']])
saveRDS(test_sce, snakemake@output[['test']])

### random subset of 500 cells
cell_ids <- sample(1:ncol(train_sce), 500)
random_subset <- train_sce[,cell_ids]

saveRDS(random_subset, snakemake@output[['random_train']])

test_labels <- tibble(cell_id = colnames(random_subset), 
                      cell_type = random_subset$CellType)

write_tsv(test_labels, 'data/scRNASeq/random/random-annotation-test.tsv')