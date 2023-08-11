suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(1)

sce <- readRDS(snakemake@input$sce)
ct_table <- tibble(cell_id = colnames(sce), cell_type = sce$cell_type)

# Cell types to select
cell_types <- snakemake@params$cell_types

# Create balanced dataset
bal_dataset <- lapply(cell_types, function(x){
  ct <- filter(ct_table, cell_type == x)
  ct <- ct[sample(1:nrow(ct), 100), ]
  ct
}) |> bind_rows()


# Create imbalanced dataset
maj_cell_type <- snakemake@params$maj_cell_type
imbalanced_dataset <- lapply(cell_types, function(x){
  # determine how many cells to select
  if(x == maj_cell_type){ # 400 if majority
    num_cells <- 400
  }else{ # 25 if minority
    num_cells <- 25
  }
  
  # select cells
  ct <- filter(ct_table, cell_type == x)
  ct <- ct[sample(1:nrow(ct), num_cells), ]
  ct
}) |> bind_rows()


minority_sce <- sce[, bal_dataset$cell_id]
majority_sce <- sce[, imbalanced_dataset$cell_id]

saveRDS(minority_sce, snakemake@output$minority)
saveRDS(majority_sce, snakemake@output$majority)

