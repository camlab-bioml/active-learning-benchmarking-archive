suppressPackageStartupMessages({
  library(tidyverse)
  library(yardstick)
})
source("pipeline/whatsthatcell-helpers.R")

log <- file(snakemake@log[[1]], open="wt")
sink(log)

sce <- readRDS(snakemake@input[['sce']])
sce_gt <- tibble(cell_id = colnames(sce),
                 annotated_cell_type = sce$CellType)

AL_entropy <- lapply(snakemake@input[['active_learning_entropy']], read_tsv) %>% 
  bind_rows()
AL_maxp <- lapply(snakemake@input[['active_learning_maxp']], read_tsv) %>%
  bind_rows()
AL <- bind_rows(AL_entropy, AL_maxp)

Seurat <- lapply(snakemake@input[['seurat']], read_tsv) %>% 
  bind_rows()
random <- lapply(snakemake@input[['random']], read_tsv) %>% 
  bind_rows()

predictions <- bind_rows(AL, Seurat, random) %>% 
  left_join(sce_gt) %>% 
  separate(prediction_params, c('m1', 'm2', 'rm_it', 'set', 'rm_knn', 'knn', 
                     'rm_res', 'res', 'rm_cell_num', 'cell_num', 'rm_rand', 'rand', 
                     'rm_corr', 'corrupted'), sep = '-') %>% 
  select(-starts_with('rm')) %>% 
  unite(method, c(m1, m2), sep = '-') %>% 
  mutate(selection_procedure = gsub("Active-Learning_entropy-strategy-", "", selection_procedure),
         selection_procedure = gsub("Active-Learning_maxp-strategy-", "", selection_procedure))

acc <- predictions %>% 
  filter(predicted_cell_type != "unassigned") %>%
  group_by(method, set, knn, res, cell_num, rand, corrupted, selection_procedure) %>% 
  acc_wrap() %>% 
  mutate(selection_procedure = gsub("_quant_entropy", "-entropy-AL", selection_procedure),
         selection_procedure = gsub("_entropy", "-entropy-AL", selection_procedure),
         selection_procedure = gsub("_quant_maxp", "-maxp-AL", selection_procedure),
         selection_procedure = gsub("_maxp", "-maxp-AL", selection_procedure),
         selection_procedure = gsub("-strategy-NA", "", selection_procedure),
         cell_num = factor(cell_num, levels = c("50", "100", "250", "500")))

write_tsv(acc, snakemake@output[['acc']])