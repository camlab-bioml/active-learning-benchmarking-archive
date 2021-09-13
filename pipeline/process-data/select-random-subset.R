library(tidyverse)
library(SingleCellExperiment)

sce <- readRDS(snakemake@input[['sce']])
sce <- readRDS("data/scRNASeq/scRNASeq-train.rds")

set.seed(42)
subset1 <- sce[, sample(1:ncol(sce), 500)]
subset2 <- sce[, sample(1:ncol(sce), 500)]
subset3 <- sce[, sample(1:ncol(sce), 500)]


gt_1 <- tibble(cell_id = colnames(subset1),
               gt_cell_type = subset1$CellType,
               method = "random",
               set = "Set1")

gt_2 <- tibble(cell_id = colnames(subset2),
               gt_cell_type = subset2$CellType,
               method = "random",
               set = "Set2")

gt_3 <- tibble(cell_id = colnames(subset3),
               gt_cell_type = subset3$CellType,
               method = "random",
               set = "Set3")


write_rds(subset1, snakemake@output[['sce_subset1']])
write_rds(subset2, snakemake@output[['sce_subset2']])
write_rds(subset3, snakemake@output[['sce_subset3']])

write_tsv(gt_1, snakemake@output[['gt_subset1']])
write_tsv(gt_2, snakemake@output[['gt_subset2']])
write_tsv(gt_3, snakemake@output[['gt_subset3']])