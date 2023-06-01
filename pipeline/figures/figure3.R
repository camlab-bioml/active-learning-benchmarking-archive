suppressPackageStartupMessages({
  library(scater)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(magick)
})
source("pipeline/whatsthatcell-helpers.R")

### [ ACCURACIES ] ####
f <- list.files("output/v8/results/", pattern = "overall", full.names = TRUE)

acc <- lapply(f, function(x){
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

schematic <- image_read("illustrator-figures/AR-schematic.ai") |> 
  image_ggplot()

eval <- full_acc_plot_wrapper(acc, "rf", "ranking", "") &
  labs(fill = "Selection method")

pdf("output/v8/paper-figures/figure3-overall-benchmark.pdf", height = 8.4, width = 9)
  eval
dev.off()

space_eval <- (plot_spacer() | eval) + 
  plot_layout(widths = c(0.001, 1))

pdf("output/v8/paper-figures/figure3.pdf", height = 12, width = 9)
  schematic / wrap_elements(full = space_eval) + 
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1.7, 4))
dev.off()




sel_acc <- filter(acc, corrupted == 0) |> 
  filter(AL_alg == "rf" | is.na(AL_alg)) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(.metric == "f_meas") |> 
  filter(initial == "ranking" | is.na(initial))
  

cross_cohort_ranks <- get_ranked_mean_estimate_by_cohort(sel_acc)
method_order <- get_cross_cohort_mean_ranked_estimate(cross_cohort_ranks) |> 
  pull(selection_procedure)

sel_acc |> 
  group_by(selection_procedure, cell_num, cohort) |> 
  mutate(selection_procedure = factor(selection_procedure, levels = method_order),
         mean_estimate = mean(.estimate),
         lower_quant = quantile(na.omit(.estimate), 0.25),
         upper_quant = quantile(na.omit(.estimate), 0.75)) |> 
  ggplot() +
  aes(fill = selection_procedure) +
  geom_segment(aes(x = as.character(cell_num),
                   y = lower_quant, yend = upper_quant),
               position = position_dodge(width = 0.9)) +
  geom_point(aes(x = as.character(cell_num), y = mean_estimate, colour = selection_procedure),
             position = position_dodge(width = 0.9)) +
  scale_color_manual(values = sel_meth_cols) +
  whatsthatcell_theme() +
  facet_wrap(~cohort, scales = "free_y", nrow = 1) +
  labs(title = metric) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


devtools::load_all("../ggplot2/")
p <- sel_acc |> 
  filter(cohort == "scRNASeq" & cell_num == 100) |> 
  group_by(selection_procedure, cell_num, cohort) |> 
  mutate(selection_procedure = factor(selection_procedure, levels = method_order),
         mean_estimate = mean(.estimate),
         lower_quant = quantile(na.omit(.estimate), 0.25),
         upper_quant = quantile(na.omit(.estimate), 0.75)) |> 
  ggplot(aes(cell_num, mean_estimate, colour = selection_procedure))

# Want the ends stay stable? Define them.
p + geom_segment(
  aes(xend = cell_num, yend = 0),
  position = position_dodge(width = 0.9)
)

p + geom_segment(aes(yend = 0), position = position_dodge(width = 0.8))
