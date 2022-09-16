library(tidyverse)
library(SingleCellExperiment)

sce <- readRDS(snakemake@input[['sce']])

set.seed(1)
subset1 <- sce[, sample(1:ncol(sce), as.integer(snakemake@wildcards[['cell_num']]))]

gt_1 <- tibble(cell_id = colnames(subset1),
               cell_type = subset1$CellType,
               method = "random",
               set = paste0("Set", snakemake@wildcards[['set_num']]),
               params = paste0("cell_num-", snakemake@wildcards[['cell_num']], "-corruption_percentage-", snakemake@wildcards[['corrupt']]))

write_rds(subset1, snakemake@output[['sce_subset1']])
write_tsv(gt_1, snakemake@output[['gt_subset1']])