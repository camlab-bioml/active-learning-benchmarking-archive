suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(magick)
  library(ggpubr)
})
source("pipeline/whatsthatcell-helpers.R")

## Part A: Proportion of cells initially selected
get_proportion_rep <- function(sce_path, selected_cells_path, cohort){
  sce <- readRDS(sce_path)
  
  cell_types <- unique(sce$cell_type)
  
  if(is.null(cell_types)){
    cell_types <- unique(sce$CellType)
  }
  
  rand_files <- list.files(selected_cells_path, 
                           pattern = "random",
                           full.names = TRUE)
  rank_files <- list.files(selected_cells_path, 
                           pattern = "ranking",
                           full.names = TRUE)
  
  random <- lapply(rand_files, function(x){
    sel_cells <- read_tsv(x) |>
      filter(iteration == 0) |> 
      pull(cell_type)
    
    sum(cell_types %in% sel_cells) / length(cell_types)
  }) |> unlist()
  
  ranking <- lapply(rank_files, function(x){
    sel_cells <- read_tsv(x) |> 
      filter(iteration == 0) |> 
      pull(cell_type)
    
    sum(cell_types %in% sel_cells) / length(cell_types)
  }) |> unlist()
  
  tibble(prop_sel = c(random, ranking), 
         type = c(rep("Random", length(random)),
                  rep("Ranked", length(ranking))),
         cohort = cohort)
}

props <- bind_rows(
  get_proportion_rep("data/CyTOF/CyTOF-train-seed-0.rds", 
                     "output/v8/data/CyTOF/Active-Learning_entropy/AL-batches-subset", 
                     "CyTOF"),
  get_proportion_rep("data/scRNASeq/scRNASeq-train-seed-0.rds", 
                     "output/v8/data/scRNASeq/Active-Learning_entropy/AL-batches-subset", 
                     "scRNASeq"),
  get_proportion_rep("data/snRNASeq/snRNASeq-train-seed-0.rds", 
                     "output/v8/data/snRNASeq/Active-Learning_entropy/AL-batches-subset", 
                     "snRNASeq")
)

prop <- props |> 
  mutate(cohort = factor(cohort, levels = c("scRNASeq", "snRNASeq", "CyTOF"))) |> 
  ggplot(aes(x = prop_sel, y = type, fill = type)) +
  geom_boxplot() +
  labs(x = "Proportion of cell types selected in initial training set of 20 cells",
       y = " \nInitial selection\nmethod") +
  scale_fill_manual(values = c("#8F3985", "#98DFEA")) +
  facet_wrap(~cohort, nrow = 1) +
  whatsthatcell_theme() + 
  theme(legend.position = "none")


## Part B: heatmaps
f <- c("output/v8/results/overall-CyTOF-benchmarking-accuracies.tsv", 
       "output/v8/results/overall-scRNASeq-benchmarking-accuracies.tsv", 
       "output/v8/results/overall-snRNASeq-benchmarking-accuracies.tsv")

acc <- lapply(f, function(x){
  df <- read_tsv(x) |> 
    mutate(cohort = case_when(grepl("CyTOF", basename(x)) ~ "CyTOF",
                              grepl("snRNASeq", basename(x)) ~ "snRNASeq",
                              grepl("scRNASeq", basename(x)) ~ "scRNASeq"))
  df
}) |> bind_rows() |>
  mutate(selection_procedure = case_when(selection_procedure == "NoMarkerSeurat_clustering" ~ "AR NoMarker",
                                         selection_procedure == "MarkerSeurat_clustering" ~ "AR Marker",
                                         selection_procedure == "0.95-entropy-AL" ~ "AL 0.95-entropy",
                                         selection_procedure == "highest-entropy-AL" ~ "AL highest-entropy",
                                         selection_procedure == "lowest-maxp-AL" ~ "AL lowest-maxp",
                                         selection_procedure == "0.05-maxp-AL" ~ "AL 0.05-maxp",
                                         selection_procedure == "0.25-maxp-AL" ~ "AL 0.25-maxp",
                                         selection_procedure == "0.75-entropy-AL" ~ "AL 0.75-entropy",
                                         TRUE ~ selection_procedure))


create_heatmap <- function(acc, sel_cohort, comp, legend = FALSE){
  subset_acc <- filter(acc, corrupted == 0) |> 
    filter(rand == 0 | is.na(rand)) |> 
    filter(.metric == 'f_meas' & cohort == sel_cohort) |> 
    filter(selection_procedure != "random" & !grepl("-AR", selection_procedure))
  
  if(comp == "lr_vs_rf"){
    improvement <- subset_acc |> 
      pivot_wider(names_from = "AL_alg", values_from = ".estimate") |> 
      mutate(improvement = (rf - multinom) / multinom) |> 
      filter(improvement != Inf & !is.na(improvement)) |> 
      filter(initial == "random") |> 
      select(-c(multinom, rf, .estimator, .metric, cohort, rand, corrupted, knn, res, initial))
  }else if(comp == "random_vs_ranked"){
    improvement <- subset_acc |> 
      pivot_wider(names_from = "initial", values_from = ".estimate") |> 
      mutate(improvement = (ranking - random) / random) |> 
      filter(improvement != Inf & !is.na(improvement)) |> 
      select(-c(ranking, random, .estimator, .metric, cohort, rand, corrupted, knn, res, AL_alg))
  }
  
  mat <- improvement |> 
    group_by(cell_num, selection_procedure) |> 
    summarize(mean_improvement = mean(improvement)) |> 
    ungroup() |> 
    mutate(mean_improvement = case_when(mean_improvement > 0.2 ~ 0.2,
                                        mean_improvement < -0.2 ~ -0.2,
                                        TRUE ~ mean_improvement)) |> 
    pivot_wider(names_from = selection_procedure, values_from = mean_improvement)
  
  ha <- HeatmapAnnotation(`Cell number` = gsub("_[0-9]", "", mat$cell_num),
                          which = "row",
                          show_legend = legend,
                          col = list(`Cell number` = c(`100` = "#FF99C9", 
                                                       `250` = "#A2C7E5", 
                                                       `500` = "#C1BDDB")))
  
  mat |> 
    as.data.frame() |> 
    column_to_rownames("cell_num") |> 
    as.matrix() |> 
    Heatmap(right_annotation = ha,
            column_title = sel_cohort,
            show_heatmap_legend = legend,
            row_order = c("100", "250", "500"),
            col = colorRamp2(seq(-0.2, 0.2, length = 3), c("blue", "#EEEEEE", "red")),
            name = "Improvement\nscore")
}

ranking_vs_random <- create_heatmap(acc, "scRNASeq", "random_vs_ranked") + 
  create_heatmap(acc, "snRNASeq", "random_vs_ranked") +
  create_heatmap(acc, "CyTOF", "random_vs_ranked", TRUE)

ranking <- draw(ranking_vs_random, column_title = "Selecting initial cells based on marker expression")

pdf("output/v8/paper-figures/ranking-random-heatmap.pdf", height = 3.5, width = 7)
  ranking
dev.off()

lr_vs_rf <- create_heatmap(acc, "scRNASeq", "lr_vs_rf") + 
  create_heatmap(acc, "snRNASeq", "lr_vs_rf") +
  create_heatmap(acc, "CyTOF", "lr_vs_rf", TRUE)

lr_vs_rf <- draw(lr_vs_rf, column_title = "F1-score improvement by random forest compared to logistic regression")

pdf("output/v8/paper-figures/lr-rf-heatmap.pdf", height = 3.5, width = 7)
  lr_vs_rf
dev.off()

## Full main figure
ranking <- image_read_pdf("output/v8/paper-figures/ranking-random-heatmap.pdf")
ranking <- ggplot() +
  background_image(ranking)

lr_vs_rf <- image_read_pdf("output/v8/paper-figures/lr-rf-heatmap.pdf")
lr_vs_rf <- ggplot() +
  background_image(lr_vs_rf)

pdf("output/v8/paper-figures/lr-and-rf-heatmaps.pdf", height = 8.8, width = 7)
   (lr_vs_rf / wrap_elements(full = prop) / ranking) + plot_layout(heights = c(3, 1.2, 3)) + plot_annotation(tag_levels = "A")
dev.off()

### Supplementary figures
# LR vs RF
rf_vs_rf_random <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(initial == "random") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  mutate(AL_alg = case_when(AL_alg == "rf" ~ "Random Forest",
                            AL_alg == "multinom" ~ "Logistic Regression")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = AL_alg)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Active learning\nalgorithm", 
       title = "Performance of selection methods comparing random forest and logistic regression active learning strategies",
       subtitle = "Shown are results obtained when randomly selecting the initial set of cells") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

rf_vs_rf_ranking <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(initial == "ranking") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  mutate(AL_alg = case_when(AL_alg == "rf" ~ "Random Forest",
                            AL_alg == "multinom" ~ "Logistic Regression")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = AL_alg)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Active learning\nalgorithm", 
       title = "Performance of selection methods comparing random forest and logistic regression active learning strategies",
       subtitle = "Shown are results obtained when ranking the initial set of cells") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))

pdf("output/v8/paper-figures/rf-vs-lr-supplementary.pdf", height = 13, width = 12)
  rf_vs_rf_random / rf_vs_rf_ranking + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
dev.off()



# random vs ranking
rand_vs_rank_lr <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(AL_alg == "multinom") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  mutate(initial = case_when(initial == "random" ~ "Random",
                            initial == "ranking" ~ "Ranking")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = initial)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Initial selection", 
       title = "Performance of selection methods comparing random and ranking based initial cell selections",
       subtitle = "Shown are results obtained using logistic regression as an active learning algorithm") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

rand_vs_rank_rf <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(AL_alg == "rf") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  mutate(initial = case_when(initial == "random" ~ "Random",
                            initial == "ranking" ~ "Ranking")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = initial)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Initial selection", 
       title = "Performance of selection methods comparing random and ranking based initial cell selections",
       subtitle = "Shown are results obtained using random forest as an active learning algorithm") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))

pdf("output/v8/paper-figures/random-vs-ranking-supplementary.pdf", height = 13, width = 12)
  rand_vs_rank_lr / rand_vs_rank_rf + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
dev.off()
