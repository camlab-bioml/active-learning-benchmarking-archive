suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
})
set.seed(1)

sce <- readRDS(snakemake@input$sce)
markers <- read_yaml(snakemake@input$markers)
percent_removed <- snakemake@wildcards$rem_percentage

unique_markers <- unique(unlist(markers$cell_types))

# Find most expressed genes
highly_expressed <- assays(sce)$logcounts |> 
  as.matrix() |> 
  rowSums() |> 
  as.data.frame() |> 
  rownames_to_column("gene") |> 
  filter(!gene %in% unique_markers)
colnames(highly_expressed)[2] <- "sum_logcounts"

highly_expressed <- arrange(highly_expressed, -sum_logcounts) |> 
  slice_head(n = 10000) |> 
  pull(gene)

# Replace varying number of markers
num_markers <- unlist(markers) |> unique() |> length()
num_markers_to_replace <- round(percent_removed * num_markers)

random_genes <- sample(highly_expressed, num_markers_to_replace)
names(random_genes) <- sample(unique(unlist(markers$cell_types)), num_markers_to_replace)


for(i in names(markers$cell_types)){
  match_idx <- which(markers$cell_types[[i]]$positive %in% names(random_genes))
  for(j in match_idx){
    markers$cell_types[[i]]$positive[j] <- random_genes[[markers$cell_types[[i]]$positive[j]]]
  }
  
  if(!is.null(markers$cell_types[[i]]$negative)){
    match_idx <- which(markers$cell_types[[i]]$negative %in% names(random_genes))
    for(j in match_idx){
      markers$cell_types[[i]]$negative[j] <- random_genes[[markers$cell_types[[i]]$negative[j]]]
    }
  }
}


write_yaml(markers, snakemake@output$markers)

