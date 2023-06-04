suppressPackageStartupMessages({
  library(scater)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(magick)
})
devtools::load_all("/ggplot2")
source("pipeline/whatsthatcell-helpers.R")

### [ ACCURACIES ] ####
acc <- lapply(snakemake@input$accs, function(x){
  df <- read_tsv(x) |> 
    mutate(cohort = case_when(grepl("CyTOF", basename(x)) ~ "CyTOF",
                              grepl("snRNASeq", basename(x)) ~ "snRNASeq",
                              grepl("scRNASeq", basename(x)) ~ "scRNASeq"))
  df
}) |> 
  bind_rows() |> 
  mutate(cohort = factor(cohort, levels = c("scRNASeq", "snRNASeq", "CyTOF")),
         selection_procedure = case_when(selection_procedure == "random" ~ "Random",
                                         selection_procedure == "NoMarkerSeurat_clustering" ~ "AR No Marker",
                                         selection_procedure == "MarkerSeurat_clustering" ~ "AR Marker",
                                         selection_procedure == "highest-entropy-AL" ~ "AL Highest-entropy",
                                         selection_procedure == "lowest-maxp-AL" ~ "AL Lowest-maxp",
                                         selection_procedure == "0.95-entropy-AL" ~ "AL 0.95-entropy",
                                         selection_procedure == "0.05-maxp-AL" ~ "AL 0.05-maxp",
                                         selection_procedure == "0.25-maxp-AL" ~ "AL 0.25-maxp",
                                         selection_procedure == "0.75-entropy-AL" ~ "AL 0.75-entropy",
                                         TRUE ~ selection_procedure))


sel_meth_cols <- sel_met_cols

eval <- full_acc_plot_wrapper(acc, "rf", "ranking", "") &
  labs(fill = "Selection method")

pdf(snakemake@output$overall_fig, height = 8.4, width = 9)
  eval
dev.off()

