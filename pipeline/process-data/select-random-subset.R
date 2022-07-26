library(tidyverse)
library(SingleCellExperiment)

sce <- readRDS(snakemake@input[['sce']])

set.seed(as.integer(snakemake@wildcards[['set_num']]))
subset1 <- sce[, sample(1:ncol(sce), as.integer(snakemake@wildcards[['cell_num']]))]
# subset2 <- sce[, sample(1:ncol(sce), as.integer(snakemake@wildcards[['cell_num']]))]
# subset3 <- sce[, sample(1:ncol(sce), as.integer(snakemake@wildcards[['cell_num']]))]

if(as.numeric(snakemake@wildcards[['corrupt']]) != 0){
    num_corrupt <- round(as.integer(snakemake@wildcards[['cell_num']]) * as.numeric(snakemake@wildcards[['corrupt']]), 0)
    all_cell_types <- unique(sce$CellType)

    corrupt_idx <- sample(1:as.integer(snakemake@wildcards[['cell_num']]), num_corrupt)

    for(i in 1:length(all_cell_types)){
        tc <- all_cell_types[i]

        to_corrupt <- which(subset1$CellType == tc)
        to_corrupt <- intersect(to_corrupt, corrupt_idx)
        if(length(to_corrupt) > 0){
            subset_cell_types <- all_cell_types[!grepl(tc, all_cell_types)]
            subset1$CellType[to_corrupt] <- sample(subset_cell_types, length(to_corrupt), replace = TRUE)
        }
    }

    # df_expression$gt_cell_type[corrupt_idx] <- sample(df_expression$gt_cell_type[corrupt_idx])

    # corrupt_idx <- sample(1:as.integer(snakemake@wildcards[['cell_num']]), num_corrupt)

    # subset1[,corrupt_idx]$CellType <- sample(subset1[,corrupt_idx]$CellType)
    # subset2[,corrupt_idx]$CellType <- sample(subset2[,corrupt_idx]$CellType)
    # subset3[,corrupt_idx]$CellType <- sample(subset3[,corrupt_idx]$CellType)
}

gt_1 <- tibble(cell_id = colnames(subset1),
               cell_type = subset1$CellType,
               method = "random",
               set = paste0("Set", snakemake@wildcards[['set_num']]),
               params = paste0("cell_num-", snakemake@wildcards[['cell_num']], "-corruption_percentage-", snakemake@wildcards[['corrupt']]))

# gt_2 <- tibble(cell_id = colnames(subset2),
#                cell_type = subset2$CellType,
#                method = "random",
#                set = "Set2",
#                params = paste0("cell_num-", snakemake@wildcards[['cell_num']], "-corruption_percentage-", snakemake@wildcards[['corrupt']]))

# gt_3 <- tibble(cell_id = colnames(subset3),
#                cell_type = subset3$CellType,
#                method = "random",
#                set = "Set3",
#                params = paste0("cell_num-", snakemake@wildcards[['cell_num']], "-corruption_percentage-", snakemake@wildcards[['corrupt']]))


write_rds(subset1, snakemake@output[['sce_subset1']])
# write_rds(subset2, snakemake@output[['sce_subset2']])
# write_rds(subset3, snakemake@output[['sce_subset3']])

write_tsv(gt_1, snakemake@output[['gt_subset1']])
# write_tsv(gt_2, snakemake@output[['gt_subset2']])
# write_tsv(gt_3, snakemake@output[['gt_subset3']])