suppressPackageStartupMessages({
  library(tidyverse)
})
source("pipeline/whatsthatcell-helpers.R")

doublets1 <- read_tsv("output/v6/results/doublet-id/scRNASeq-doublet-id-10x Chromium (v2) A.tsv")
colnames(doublets1)[2:3] <- c("score", "classification")

doublets2 <- read_tsv("output/v6/results/doublet-id/scRNASeq-doublet-id-10x Chromium (v2) B.tsv")
colnames(doublets2)[2:3] <- c("score", "classification")

doublets <- bind_rows(
  doublets1,
  doublets2
)

files <- list.files("output/v6/results/rem_cell_type/", full.names = TRUE)
entropies <- lapply(files, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('AL_alg-multinom', params) ~ "multinom",
                      grepl('AL_alg-rf', params) ~ 'rf'),
       strat = case_when(grepl('strat-highest_entropy', params) ~ 'highest_entropy',
                         grepl('strat-lowest_maxp', params) ~ 'lowest_maxp'),
       init = case_when(grepl('init-random', params) ~ "random",
                        grepl('init-ranking', params) ~ "ranking"),
       ct = str_extract(params, 'rem_celltype.*'),
       s = str_extract(params, '-seed-.*')) |>
  mutate(ct = gsub("rem_celltype-", "", ct),
         ct = gsub("-seed-[0-9]", "", ct),
         s = gsub("-seed-", "", s)) |> 
  select(-params)

pdf(snakemake@output$pdf, height = 8, width = 20)
  entropies |> 
    left_join(select(doublets, -params), by = "cell_id") |> 
    ggplot(aes(x = gt_cell_type, y = criterion_val, fill = classification)) +
    geom_boxplot() +
    labs(x = "Ground truth", y = "Entropy", "Doublet classification") +
    facet_grid(num_missing_cells ~ ct) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()



