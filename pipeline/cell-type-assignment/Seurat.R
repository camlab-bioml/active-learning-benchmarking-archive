library(Seurat)
library(SingleCellExperiment)
library(yaml)
library(tidyverse)
library(glue)

### [READ IN FILES] ###
sce <- readRDS(snakemake@input[['training_rds']])
markers <- read_yaml(snakemake@input[['markers']])$cell_types
unique_markers <- unlist(markers) %>% unique()


### [PROCESS DATA] ###
if(length(names(assays(sce))) == 2){
  seu <- CreateSeuratObject(counts = assays(sce)$counts)
} else{
  seu <- CreateSeuratObject(counts = assays(sce)$logcounts)
}
seu <- NormalizeData(seu)

seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

### [CLUSTERING] ###
seu <- FindNeighbors(seu, dims = 1:snakemake@params[['max_dim']])
seu <- FindClusters(seu, resolution  = snakemake@params[['resolution']])

# First UMAP viz
seu <- RunUMAP(seu, dims = 1:10)
pdf(snakemake@output[['cluster_umap_pdf']])
DimPlot(seu, reduction = 'umap')
dev.off()

# Get clusters
clusters <- Idents(seu) %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id') %>% 
  rename(cluster = '.')

# Get expression
expression <- seu@assays$RNA@counts %>% 
  as.matrix() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')

expression <- expression[, c('cell_id', unique_markers)] %>% 
  left_join(clusters) %>% 
  select(-cell_id)


# Calculate enrichments for each cell type/cluster
enrichments <- lapply(1:length(markers), function(x){
  p <- expression[, c(markers[[x]]$positive, 'cluster')] %>% 
    pivot_longer(-cluster, names_to = "marker", values_to = "expression") %>% 
    group_by(cluster) %>% 
    summarize(mean = mean(expression))
  
  if(!is.null(markers[[x]]$negative)){
    n <- expression[, c(markers[[x]]$negative, 'cluster')] %>% 
      pivot_longer(-cluster, names_to = "marker", values_to = "expression") %>% 
      group_by(cluster) %>% 
      summarize(mean = mean(expression))
    
    diff <- p$mean - n$mean
  }else{
    diff <- p$mean
  }
  
  assignment <- tibble(cluster = p$cluster, 
                       avrg_expression = diff,
                       cell_type = names(markers[x]))
  assignment
}) %>% bind_rows()

# Find cell type with max value for each cluster
cluster_assignments <- enrichments %>% 
  group_by(cluster) %>% 
  mutate(max_score = avrg_expression == max(avrg_expression)) %>%
  filter(max_score == TRUE)

# Join cell id's with cell type label
assignments <- left_join(clusters, cluster_assignments) %>% 
  select(cell_id, cell_type, cluster) %>% 
  rename(predicted_cell_type = cell_type) %>%
  mutate(prediction_params = paste0("Seurat - Cluster annotation - Average expression - ",
                                    "max_pca_dim ", snakemake@wildcards[['max_pca_dim']],
                                    " - resolution ", snakemake@wildcards[['res']]),
         modality = snakemake@wildcards[['modality']])

assignments %>% 
  select(cell_id, predicted_cell_type, prediction_params, cluster, modality) %>% 
  write_tsv(snakemake@output[['assignments']])

### Create diagnostic plots
# Positive markers
lapply(1:length(markers), function(x){
  cell_types <- names(markers[x]) %>% 
    gsub(" ", "-", .)
  
  out_path <- snakemake@params[['positive_markers_diagnostic']] %>% 
    gsub("\\[", "{", .) %>% 
    gsub("\\]", "}", .)
  out <- glue(out_path)
  
  pdf(out, width = 12)
  print(FeaturePlot(seu, features = markers[[x]]$positive))
  dev.off()
})

# Negative markers
lapply(1:length(markers), function(x){
  cell_types <- names(markers[x]) %>% 
    gsub(" ", "-", .)
  
  out_path <- snakemake@params[['negative_markers_diagnostic']] %>% 
    gsub("\\[", "{", .) %>% 
    gsub("\\]", "}", .)
  out <- glue(out_path)
  
  if(!is.null(markers[[x]]$negative)){
    pdf(out, width = 12)
    print(FeaturePlot(seu, features = markers[[x]]$negative))
    dev.off()
  }
})


### Cell type umap
Idents(seu) <- assignments$predicted_cell_type

pdf(snakemake@output[['cell_type_umap_pdf']])
DimPlot(seu, reduction = 'umap', label = TRUE)
dev.off()

Idents(seu) <- sce$CellType
pdf(snakemake@output[['ground_truth_umap_pdf']])
DimPlot(seu, reduction = 'umap', label = TRUE)
dev.off()