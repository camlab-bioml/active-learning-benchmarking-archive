suppressPackageStartupMessages(
  library(tidyverse)
)
source("pipeline/whatsthatcell-helpers.R")


cytof_acc <- read_tsv("output/v7/results/overall-CyTOF-benchmarking-accuracies.tsv") |> 
  mutate(cohort = "CyTOF")
scrna_acc <- read_tsv("output/v7/results/overall-scRNASeq-benchmarking-accuracies.tsv") |> 
  mutate(cohort = "scRNASeq")
snrna_acc <- read_tsv("output/v7/results/overall-snRNASeq-benchmarking-accuracies.tsv") |> 
  mutate(cohort = "snRNASeq")

acc <- bind_rows(cytof_acc, scrna_acc, snrna_acc)



plot_knn_res <- function(acc, fill, cohort){
  acc <- filter(acc, selection_procedure == "MarkerSeurat_clustering" |
                  selection_procedure == "NoMarkerSeurat_clustering") |> 
    mutate(selection_procedure = case_when(selection_procedure == "MarkerSeurat_clustering" ~ "AR Marker",
                                           selection_procedure == "NoMarkerSeurat_clustering" ~ "AR No Marker")) |> 
    filter(.metric == "f_meas")
  
  if(fill == "res"){
    acc |>
      ggplot(aes(x = as.character(cell_num), y = .estimate, fill = as.character(res))) +
      geom_boxplot() +
      labs(x = "Number of cells", y = "F1 score", fill = "Clustering\nresolution",
           title = cohort) +
      facet_grid(method ~ selection_procedure + knn) +
      whatsthatcell_theme()
  }else if(fill == "knn"){
    acc |> 
      ggplot(aes(x = as.character(cell_num), y = .estimate, fill = as.character(knn))) +
      geom_boxplot() +
      labs(x = "Number of cells", y = "F1 score", fill = "Number of\nnearest\nneighbours",
           title = cohort) +
      facet_grid(method ~ selection_procedure + res) +
      whatsthatcell_theme()
  }
}


pdf("output/v7/paper-figures/supp-AR-res.pdf", height = 14, width = 9)
  (plot_knn_res(cytof_acc, "res", "CyTOF") /
    plot_knn_res(scrna_acc, "res", "scRNASeq") /
    plot_knn_res(snrna_acc, "res", "snRNASeq")) +
    plot_layout(guides = "collect", heights = c(1, 2, 2))
dev.off()


pdf("output/v7/paper-figures/supp-AR-knn.pdf", height = 14, width = 9)
(plot_knn_res(cytof_acc, "knn", "CyTOF") /
    plot_knn_res(scrna_acc, "knn", "scRNASeq") /
    plot_knn_res(snrna_acc, "knn", "snRNASeq")) +
  plot_layout(guides = "collect", heights = c(1, 2, 2))
dev.off()



