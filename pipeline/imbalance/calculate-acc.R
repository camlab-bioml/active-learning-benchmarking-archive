suppressPackageStartupMessages({
  library(yardstick)
  library(tidyverse)
})
source("pipeline/whatsthatcell-helpers.R")
save.image('debug-imbalance.Rdata')
balanced_different <- list.files(snakemake@params[['balanced_diff']], full.names = TRUE)
balanced_similar <- list.files(snakemake@params[['balanced_sim']], full.names = TRUE)

imbalanced_different <- list.files(snakemake@params[['imbalanced_diff']], full.names = TRUE)
imbalanced_similar <- list.files(snakemake@params[['imbalanced_sim']], full.names = TRUE)

predictions <- lapply(list(balanced_different, balanced_similar, imbalanced_different, imbalanced_similar), function(x){
    lapply(x, function(y){
        read_tsv(y)
    }) |>
    bind_rows()
}) |> bind_rows()

gt <- lapply(snakemake@input[['sces']], function(x){
    sce <- readRDS(x)
    gt <- tibble(cell_id = colnames(sce),
             annotated_cell_type = sce$CellType)

    gt
}) |> bind_rows()

predictions_gt <- left_join(predictions, gt) |>
    separate(prediction_params, c('m1', 'm2', 'rm_it', 'set', 'rm_knn', 'knn', 
                     'rm_res', 'res', 'rm_cell_num', 'cell_num', 'rm_rand', 'rand', 
                     'rm_corr', 'corrupted'), sep = '-') |>
  unite(method, c(m1, m2), sep = '-') |> 
  mutate(selection_procedure = gsub("Active-Learning_entropy-strategy-", "", selection_procedure),
         selection_procedure = gsub("Active-Learning_maxp-strategy-", "", selection_procedure),
         selection_procedure = gsub("-strategy-NA", "", selection_procedure),
         selection_procedure = gsub("Seurat-clustering", "Seurat_clustering", selection_procedure)) |> 
  separate(selection_procedure, c('strat', 'rm_AL', 'al'), sep = '-') |> 
  select(-starts_with('rm'))

acc <- predictions_gt |>
  filter(predicted_cell_type != "unassigned") |> 
  group_by(method, set, knn, res, cell_num, rand, corrupted, strat, 
           al, similarity, modality) |> 
  acc_wrap()

write_tsv(acc, snakemake@output[['acc']])
