suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(data.table)
  library(RPortfolioSimilarity)
  library(ComplexHeatmap)
})
source("pipeline/whatsthatcell-helpers.R")


sce <- readRDS("data/scRNASeq/scRNASeq-train-seed-0.rds")

df_expression <- load_scs(sce) |> 
  select(-X1)

# Remove genes with 0 expression
expression <- df_expression[, colSums(df_expression) > 0]

# Run PCA
df_PCA <- as.matrix(expression) |> 
  prcomp(center = TRUE, scale. = TRUE)

# Get variance explained
var_explained <- summary(df_PCA)$importance |> 
  t() |>
  as.data.frame() |> 
  rownames_to_column("PC") |> 
  select(PC, "Proportion of Variance") 
var_explained <- var_explained[1:20,]

# Get PCA embedding
df_PCA <- df_PCA$x[,1:20] |> 
  as.data.frame() |> 
  mutate(cell_type = sce$CellType)


celltypes <- unique(df_PCA$cell_type)

cell_type_sims <- lapply(celltypes, function(cell_type_i){
  sim <- lapply(celltypes, function(cell_type_j){
    cell_type_i_pca <- filter(df_PCA, cell_type == cell_type_i) |> 
      select(-cell_type)
    
    cell_type_j_pca <- filter(df_PCA, cell_type == cell_type_j) |> 
      select(-cell_type)
    
    cell_type_i_avg <- colSums(cell_type_i_pca) / nrow(cell_type_i_pca)
    cell_type_j_avg <- colSums(cell_type_j_pca) / nrow(cell_type_j_pca)
    
    wtVCosSimilarity(cell_type_i_avg, cell_type_j_avg, 
                     var_explained$`Proportion of Variance`)
  }) |> unlist()
  names(sim) <- celltypes
  sim
}) |> bind_rows()



cell_type_sims <- as.data.frame(cell_type_sims)
rownames(cell_type_sims) <- celltypes
Heatmap(cell_type_sims)


