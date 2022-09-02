library(caret)
library(data.table)
library(yaml)
library(SingleCellExperiment)
library(tidyverse)
source("pipeline/whatsthatcell-helpers.R")
set.seed(42)

### [ PARAMS ] #####
max_AL_iterations <- as.integer(snakemake@params[['max_cell_num']]) / 10

### [ LOAD & PROCESS DATA ] #####
markers <- read_yaml(snakemake@input[['markers']])
unique_markers <- unique(unlist(markers$cell_types))

sce <- readRDS(snakemake@input[['expression']])
if(TRUE %in% grepl("CD16_32", rownames(sce))){
  rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
}

df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$gt_cell_type <- sce$CellType
df_expression$corrupted <- NA

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
      subset_cell_types <- all_cell_types[!(tc == all_cell_types)]
      df_expression$corrupted[to_corrupt] <- sample(subset_cell_types, 
                                                    length(to_corrupt), 
                                                    replace = TRUE)
    }
  }
  
  df_expression$gt_cell_type[!is.na(df_expression$corrupted)] <- 
    df_expression$corrupted[!is.na(df_expression$corrupted)]
  df_expression$corrupted <- NULL
}

# Get list of cell types
all_cell_types <- sce$CellType |> 
  unique()

# Get number of iterations to run ranked assignment for
# iterations <- ((nrow(df_expression) - 2*length(all_cell_types)) / length(all_cell_types)) - 1

### [ RANKED CELL TYPE ASSIGNMENT ] #####
# for(i in 1:iterations){
#   # Get initial set of cells based on their marker expression ranking
#   df_expression <- cell_ranking_wrapper(df_expression, markers)
#   
#   if(all(all_cell_types %in% unique(df_expression$cell_type))){
#     print('ranking done')
#     break
#   }
# }
df_expression <- cell_ranking_wrapper(df_expression, markers)

### [ ACTIVE LEARNING CELL TYPE ASSIGNMENT ] #####
if(snakemake@wildcards[['AL_type']] == 'Active-Learning_maxp'){
  criterion <- "maxp"
}else if(snakemake@wildcards[['AL_type']] == 'Active-Learning_entropy'){
  criterion <- "entropy"
}

# create a list to save all entropies into
entropies <- list()

# Remove genes with 0 expression
df_expression <- df_expression[, c(TRUE, 
                                   colSums(df_expression[,2:(ncol(df_expression)-4)]) > 0,
                                   rep(TRUE, 4))]

# Calculate PCA embedding
df_PCA <- select(df_expression, -c(X1, cell_type, iteration, gt_cell_type, corrupted)) |> 
  as.matrix() |> 
  prcomp(center = TRUE, scale. = TRUE)

df_PCA <- df_PCA$x |> 
  as.data.frame()

df_PCA <- bind_cols(
  tibble(X1 = df_expression$X1),
  df_PCA[,1:min(20, ncol(df_PCA))], 
  tibble(cell_type = df_expression$cell_type,
         gt_cell_type = df_expression$gt_cell_type,
         iteration = NA)
)

for(i in 1:max_AL_iterations){
  AL <- active_learning_wrapper(select(df_PCA, -gt_cell_type, -iteration), 
                                snakemake@wildcards[['AL_alg']], 
                                snakemake@wildcards[['strat']], 
                                i, 
                                entropies, 
                                as.numeric(snakemake@wildcards[['rand']]),
                                criterion)

  entropies[[length(entropies) + 1]] <- AL$criterion_table
  
  # What index do the selected cells correspond to?
  to_assign_index <- match(AL$new_cells, df_PCA$X1)
  
  # Get ground truth labels based on the index
  df_PCA$cell_type[to_assign_index] <- df_PCA$gt_cell_type[to_assign_index]
  df_PCA$iteration[to_assign_index] <- i
  
  not_annotated <- filter(df_PCA, is.na(cell_type)) %>% 
    nrow()
  print(i)
  if(not_annotated < 10){
    break
  }
}


### [ SAVE OUTPUT ] #####
entropies %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(method = paste0("Active-Learning-groundTruth-strategy-", 
                         snakemake@wildcards[['strat']], '-AL_alg-', 
                         snakemake@wildcards[['AL_alg']], '-randomCells-', 
                         snakemake@wildcards[['rand']], '-corrupted-', 
                         snakemake@wildcards[['corrupt']])) %>% 
  write_tsv(snakemake@output[['entropy']])

original_cell_types <- tibble(cell_id = colnames(sce),
                              cell_type = sce$CellType)

df_PCA %>% 
  select(X1, cell_type, iteration) %>%
  dplyr::rename("cell_id" = "X1", 
                "corrupted_cell_type" = "cell_type") %>% # this is the field iteratively filled in by AL, thus may be corrupted
  left_join(original_cell_types, by = 'cell_id') %>% 
  filter(!is.na(corrupted_cell_type)) %>% 
  mutate(method = paste0("Active-Learning-groundTruth-strategy-", 
                         snakemake@wildcards[['strat']], '-AL_alg-', 
                         snakemake@wildcards[['AL_alg']], '-randomCells-', 
                         snakemake@wildcards[['rand']], '-corrupted-', 
                         snakemake@wildcards[['corrupt']])) %>% 
  write_tsv(snakemake@output[['assignments']])