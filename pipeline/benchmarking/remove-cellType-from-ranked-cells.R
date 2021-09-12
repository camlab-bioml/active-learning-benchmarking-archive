library(tidyverse)
library(yaml)
library(data.table)
library(caret)
source("pipeline/whatsthatcell-helpers.R")


markers <- read_yaml(snakemake@input[['markers']])
unique_markers <- unique(unlist(markers$cell_types))

sce <- readRDS(snakemake@input[['expression']])
df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$iteration <- NA
df_expression$gt_cell_type <- sce$CellType

### Create ground truth tibble
ground_truth <- tibble(cell_id = rownames(colData(sce)),
                       cell_type = sce$CellType)


remove_marker_AL <- function(markers, remove_index, df_expression, ground_truth){
  subset_markers <- markers$cell_types[-remove_index]
  
  # Get list of cell types
  all_cell_types <- names(subset_markers)
  
  iterations <- ((nrow(df_expression) - 2*length(all_cell_types)) / length(all_cell_types))
  iterations <- floor(iterations) - 1
  
  for(i in 1:iterations){
    # Get initial set of cells based on their marker expression ranking
    ranked_cells <- select_initial_cells(df_expression, subset_markers)
    
    # What index do the selected cells correspond to?
    to_assign_index <- match(ranked_cells, df_expression$X1)
    
    # Get ground truth labels based on the index
    assignment <- df_expression$gt_cell_type[to_assign_index]
    
    df_expression$cell_type[to_assign_index] <- assignment
    df_expression$iteration[to_assign_index] <- 0

    if(all(all_cell_types %in% unique(df_expression$cell_type))){
      break
    }
  }
  
  
  # If any cell types have been assigned to the removed cell type make them NA again
  df_expression$cell_type[which(df_expression$cell_type == names(markers$cell_types[remove_index]))] <- NA
  
  for(i in 1:20){
    # AL selected cells
    AL <- select_cells_classifier(df_expression, unique_markers)
    new_cells <- AL$selected_cells
    
    # What index do the selected cells correspond to?
    to_assign_index <- match(new_cells, df_expression$X1)
    # Get ground truth labels based on the index
    assignment <- ground_truth$cell_type[match(new_cells, ground_truth$cell_id)]
    
    df_expression$cell_type[to_assign_index] <- assignment
    df_expression$iteration[to_assign_index] <- i
  }
  
  df_expression %>% 
    filter(!is.na(cell_type)) %>% 
    select(X1, iteration, cell_type, gt_cell_type) %>% 
    mutate(removed_markers = names(markers$cell_types[remove_index]))
}


AL_removed <- lapply(1:length(markers$cell_types), function(x){
  remove_marker_AL(markers, x, df_expression, ground_truth)
}) %>% 
  bind_rows()

write_tsv(AL_removed, snakemake@output[['tsv']])

# AL_removed %>%
#   bind_rows() %>% 
#   group_by(removed_markers, iteration, cell_type) %>% 
#   tally() %>% 
#   ggplot(aes(x = iteration, y = n, fill = cell_type)) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = cell_type_colours()) +
#   facet_wrap(~removed_markers) +
#   whatsthatcell_theme()
