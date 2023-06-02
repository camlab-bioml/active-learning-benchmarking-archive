library(tidyverse)
library(ComplexHeatmap)

cytof_sim <- read_tsv("output/v7/results/cell-type-similarity/similarity-CyTOF-seed-0.tsv")
scrna_sim <- read_tsv("output/v7/results/cell-type-similarity/similarity-scRNASeq-seed-0.tsv")
snrna_sim <- read_tsv("output/v7/results/cell-type-similarity/similarity-snRNASeq-seed-0.tsv")

plot_dist_mat <- function(df){
  pivot_wider(df, names_from = "cell_type2", values_from = "cosine_similarity") |> 
    select(-c(cohort, seed)) |> 
    as.data.frame() |> 
    column_to_rownames("cell_type1") |> 
    as.matrix() |> 
    Heatmap(name = "Cosine\nsimilarity")
}

pdf("output/v7/paper-figures/Supp-cell-type-sim-scRNASeq.pdf", height = 5, width = 5.5)
  plot_dist_mat(scrna_sim)
dev.off()

pdf("output/v7/paper-figures/Supp-cell-type-sim-snRNASeq.pdf", height = 5, width = 5.5)
  plot_dist_mat(snrna_sim)
dev.off()

pdf("output/v7/paper-figures/Supp-cell-type-sim-CyTOF.pdf", height = 5, width = 5.5)
  plot_dist_mat(cytof_sim)
dev.off()
