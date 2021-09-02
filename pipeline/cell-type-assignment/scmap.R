suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scmap)
  library(tidyverse)
})

# Read in training data expression
train_sce <- readRDS(snakemake@input[['train_data']])

# Read in annotated labels
annotations <- read_tsv(snakemake@input[['annotation']]) %>% 
  as.data.frame() %>% 
  column_to_rownames('cell_id')

# Add labels to expression
train_sce$cell_type1 <- annotations[colnames(train_sce), 'cell_type']

train_sce <- train_sce[,!is.na(train_sce$cell_type1)]

# Process data
rowData(train_sce)$feature_symbol <- rownames(train_sce)
train_sce <- selectFeatures(train_sce, suppress_plot = TRUE)

# Read in test set
annotate_sce <- readRDS(snakemake@input[['test_data']])
rowData(annotate_sce)$feature_symbol <- rownames(annotate_sce)

### Start training
train_sce <- indexCluster(train_sce)


# Map clusters
scmapCluster_results <- scmapCluster(
  projection = annotate_sce, 
  index_list = list(
    yan = metadata(train_sce)$scmap_cluster_index
  )
)

clustering_prediction <- tibble(cell_id = colnames(annotate_sce),
                                predicted_cell_type = scmapCluster_results$combined_labs,
                                prediction_params = paste0('scmap-clusters-iterations_set-', 
                                                           snakemake@wildcards[['set']],
                                                           "max_dim_pca", snakemake@wildcards[['pca']],
                                                           "res", snakemake@wildcards[['res']]),
                                selection_procedure = snakemake@wildcards[['selection_procedure']],
                                training_annotator = snakemake@wildcards[['annotator']],
                                modality = 'scRNASeq')


write_tsv(clustering_prediction, snakemake@output[['cluster_predictions']])

### Cell level
set.seed(1)

train_sce <- indexCell(train_sce)

scmapCell_results <- scmapCell(
  annotate_sce, 
  list(
    yan = metadata(train_sce)$scmap_cell_index
  )
)

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(train_sce)$cell_type1)
  )
)

sc_prediction <- tibble(cell_id = colnames(annotate_sce),
                        predicted_cell_type = scmapCell_clusters$scmap_cluster_labs[,1],
                        prediction_params = paste0('scmap-sc-iterations_set-',
                                                   snakemake@wildcards[['set']],
                                                   "-max_dim_pca-", snakemake@wildcards[['pca']],
                                                   "-res-", snakemake@wildcards[['res']]),
                        selection_procedure = snakemake@wildcards[['selection_procedure']],
                        training_annotator = snakemake@wildcards[['annotator']],
                        modality = 'scRNASeq')

write_tsv(sc_prediction, snakemake@output[['sc_predictions']])

