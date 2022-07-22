library(caret)
library(data.table)
library(yaml)
library(SingleCellExperiment)
library(tidyverse)
source("pipeline/whatsthatcell-helpers.R")
#save.image('debug-AL')
set.seed(42)

### [ LOAD & PROCESS DATA ] #####
markers <- read_yaml(snakemake@input[['markers']])
unique_markers <- unique(unlist(markers$cell_types))

sce <- readRDS(snakemake@input[['expression']])
if(TRUE %in% grepl("CD16_32", rownames(sce))){
  rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
}

df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$iteration <- NA
df_expression$gt_cell_type <- sce$CellType

## Corrupt a fraction of labels
if(as.numeric(snakemake@wildcards[['corrupt']]) != 0){
  all_cell_types <- unique(sce$CellType)
  num_corrupt <- round((nrow(df_expression) * as.numeric(snakemake@wildcards[['corrupt']])), 0)
  corrupt_idx <- sample(1:nrow(df_expression), num_corrupt)

  for(i in 1:length(all_cell_types)){
    tc <- all_cell_types[i]

    to_corrupt <- which(df_expression$gt_cell_type == tc)
    to_corrupt <- intersect(to_corrupt, corrupt_idx)
    if(length(to_corrupt) > 0){
      subset_cell_types <- all_cell_types[!grepl(tc, all_cell_types)]
      df_expression$gt_cell_type[to_corrupt] <- sample(subset_cell_types, 
                                                       length(to_corrupt), 
                                                       replace = TRUE)
    }
  }

  # df_expression$gt_cell_type[corrupt_idx] <- sample(df_expression$gt_cell_type[corrupt_idx])
}

# Get list of cell types
all_cell_types <- df_expression$gt_cell_type %>% unique()

# Get number of iterations to run ranked assignment for
iterations <- ((nrow(df_expression) - 2*length(all_cell_types)) / length(all_cell_types)) - 1


### [ RANKED CELL TYPE ASSIGNMENT ] #####
for(i in 1:iterations){
  # Get initial set of cells based on their marker expression ranking
  df_expression <- cell_ranking_wrapper(df_expression, markers)
  
  if(all(all_cell_types %in% unique(df_expression$cell_type))){
    print('ranking done')
    break
  }
}

### [ ACTIVE LEARNING CELL TYPE ASSIGNMENT ] #####

# create a list to save all entropies into
entropies <- list()
for(i in 1:nrow(df_expression)){
  AL <- active_learning_wrapper(df_expression, unique_markers, snakemake@wildcards[['strat']], i, entropies, as.numeric(snakemake@wildcards[['rand']]))

  df_expression <- AL$expression
  entropies <- AL$entropies
  
  not_annotated <- filter(df_expression, is.na(cell_type)) %>% 
    nrow()
  if(not_annotated < 10){
    break
  }
}


### [ SAVE OUTPUT ] #####
entropies %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  write_tsv(snakemake@output[['entropy']])

df_expression %>% 
  filter(!is.na(cell_type)) %>% 
  dplyr::rename("cell_id" = "X1",
                "corrupted_AL_cell_type" = "cell_type",
                "cell_type" = "gt_cell_type") %>% 
  mutate(method = paste0("Active-Learning-groundTruth-strategy-", snakemake@wildcards[['strat']], 
                         '-randomCells-', snakemake@wildcards[['rand']], '-corrupted-', snakemake@wildcards[['corrupt']])) %>% 
  select(cell_id, cell_type, method, iteration) %>% 
  write_tsv(snakemake@output[['assignments']])

