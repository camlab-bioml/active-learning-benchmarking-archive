suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

set.seed(as.integer(snakemake@wildcards$s))

sce <- readRDS(snakemake@input[['sce']])
save.image('debug-imbalance.Rdata')

majority_ct <- snakemake@params[['majority']]
minority_ct <- snakemake@params[['minority']]

if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}


# subset to majority and minority cell type datasets & sample
big_sce <- sce[, sce$CellType == majority_ct]
big_sce <- big_sce[, sample(1:ncol(big_sce), 1000)]

small_sce <- sce[, sce$CellType == minority_ct]
small_sce <- small_sce[, sample(1:ncol(small_sce), 100)]

comb_sce <- cbind(big_sce, small_sce)

# Create test set
rem_sce <- sce[, !(colnames(sce) %in% colnames(comb_sce))]
rem_sce <- rem_sce[, rem_sce$CellType == majority_ct | rem_sce$CellType == minority_ct]

saveRDS(comb_sce, snakemake@output[['sce']])
saveRDS(rem_sce, snakemake@output[['rem_sce']])