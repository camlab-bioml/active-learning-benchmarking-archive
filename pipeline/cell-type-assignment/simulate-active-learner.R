library(caret)
library(data.table)
library(yaml)
library(SingleCellExperiment)
library(tidyverse)
source("pipeline/whatsthatcell-helpers.R")


### [ LOAD & PROCESS DATA ] #####
markers <- read_yaml('markers/scRNASeq.yml')#snakemake@input[['markers']])
unique_markers <- unique(unlist(markers$cell_types))

sce <- readRDS("data/scRNASeq/scRNASeq-train.rds")#snakemake@input[['expression']])
df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$iteration <- NA
df_expression$gt_cell_type <- sce$CellType

### Create ground truth tibble
#ground_truth <- tibble(cell_id = rownames(colData(sce)),
#                       cell_type = sce$CellType)

# Get list of cell types
all_cell_types <- df_expression$gt_cell_type %>% unique()

# Get number of iterations to run ranked assignment for
iterations <- ((nrow(df_expression) - 2*length(all_cell_types)) / length(all_cell_types)) - 1


### [ RANKED CELL TYPE ASSIGNMENT ] #####
for(i in 1:iterations){
  # Get initial set of cells based on their marker expression ranking
  ranked_cells <- select_initial_cells(df_expression, markers$cell_types)
  
  # What index do the selected cells correspond to?
  to_assign_index <- match(ranked_cells, df_expression$X1)
  
  df_expression$cell_type[to_assign_index] <- df_expression$gt_cell_type[to_assign_index]
  # Get ground truth labels based on the index
  #assignment <- ground_truth$cell_type[match(ranked_cells, ground_truth$cell_id)]
  #df_expression$cell_type[to_assign_index] <- assignment
  df_expression$iteration <- 0
  
  if(all(all_cell_types %in% unique(df_expression$cell_type))){
    break
  }
}

### [ ACTIVE LEARNING CELL TYPE ASSIGNMENT ] #####

# create a list to save all entropies into
entropies <- list()
for(i in 1:nrow(df_expression)){
  # AL selected cells
  AL <- select_cells_classifier(df_expression, unique_markers)
  new_cells <- AL$selected_cells
  entropies[[length(entropies) + 1]] <- AL$entropy_table
  
  # What index do the selected cells correspond to?
  to_assign_index <- match(new_cells, df_expression$X1)
  # Get ground truth labels based on the index
  df_expression$cell_type[to_assign_index] <- df_expression$gt_cell_type[to_assign_index]
  #assignment <- ground_truth$cell_type[match(new_cells, ground_truth$cell_id)]
  
  #df_expression$cell_type[to_assign_index] <- assignment
  df_expression$iteration[to_assign_index] <- i
  
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
  dplyr::rename("cell_id" = "X1") %>% 
  mutate(method = "Active-Learning - ground truth") %>% 
  select(cell_id, cell_type, method, iteration) %>% 
  write_tsv(snakemake@output[['assignments']])

