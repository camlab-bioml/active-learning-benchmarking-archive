suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("pipeline/whatsthatcell-helpers.R")

# Read in data
### scRNASeq
scrna <- c(
  list.files("output/v7/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-JIMT1*highest*multi*20*"),
             full.names = TRUE),
  list.files("output/v7/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-AU565*highest*multi*20*"),
             full.names = TRUE),
  list.files("output/v7/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-JIMT1*highest*rf*20*"),
             full.names = TRUE),
  list.files("output/v7/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-AU565*highest*rf*20*"),
             full.names = TRUE)
)

scrna <- lapply(scrna, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))


scrna_box <- scrna |> 
  mutate(al = case_when(al == "multinom" ~ "LR",
                        al == "rf" ~ "RF")) |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_grid(al ~ params) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


### snRNASeq
snrna <- c(
  list.files("output/v7/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-Ductal*highest*multi*20*"),
             full.names = TRUE),
  list.files("output/v7/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-Endothelial*highest*multi*20*"),
             full.names = TRUE),
  list.files("output/v7/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-Schwann*highest*multi*20*"),
             full.names = TRUE),
  list.files("output/v7/results/rem_cell_type_group/",
             pattern = glob2rx("*random*l1*sn*high*multi*num-20*"),
             full.names = TRUE)
)

snrna <- lapply(snrna, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))

snrna_box <- snrna |> 
  mutate(params = factor(params, levels = c("Ductal", "Endothelial", "Schwann", 
                                            "Endothelial, Schwann"))) |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

### Cytof
cytof <- c(
  list.files("output/v7/results/rem_cell_type",
             pattern = glob2rx("*random*rem-Classical*CyTOF*highest*multi*20*"),
             full.names = TRUE),
  list.files("output/v7/results/rem_cell_type",
             pattern = glob2rx("*random*CD8*CyTOF*highest*multi*20*"),
             full.names = TRUE)
)


cytof <- lapply(cytof, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))

cytof_box <- cytof |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params, ncol = 1) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


## Combine
pdf("output/v8/paper-figures/rem-cell-type.pdf", height = 6, width = 14)
(scrna_box | snrna_box | cytof_box) +
  plot_layout(guides = "collect", widths = c(2, 2, 1)) +
  plot_annotation(tag_levels = "A") & 
  theme(legend.position = 'bottom', plot.tag = element_text(size = 18))
dev.off()



## Supplemental for RF for snRNASeq and CyTOF
snrna_sup <- c(
  list.files("output/v8/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-Ductal*highest*rf*20*"),
             full.names = TRUE),
  list.files("output/v8/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-Endothelial*highest*rf*20*"),
             full.names = TRUE),
  list.files("output/v8/results/rem_cell_type/", 
             pattern = glob2rx("Init-sel-random-rem-Schwann*highest*rf*20*"),
             full.names = TRUE),
  list.files("output/v8/results/rem_cell_type_group/",
             pattern = glob2rx("*random*l1*sn*high*rf*num-20*"),
             full.names = TRUE)
)

snrna_sup <- lapply(snrna_sup, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))

cytof_sup <- c(
  list.files("output/v8/results/rem_cell_type",
             pattern = glob2rx("*random*rem-Classical*CyTOF*highest*rf*20*"),
             full.names = TRUE),
  list.files("output/v8/results/rem_cell_type",
             pattern = glob2rx("*random*CD8*CyTOF*highest*rf*20*"),
             full.names = TRUE)
)


cytof_sup <- lapply(cytof_sup, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))


snrna_sup_box <- snrna_sup |> 
  mutate(params = factor(params, levels = c("Ductal", "Endothelial", "Schwann", 
                                            "Endothelial, Schwann"))) |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



cytof_sup_box <- cytof_sup |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params, ncol = 1) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf("output/v8/paper-figures/supp-rem-cell-type.pdf", height = 6, width = 10)
  (snrna_sup_box | cytof_sup_box) +
    plot_layout(guides = "collect", widths = c(2, 1)) +
    plot_annotation(tag_levels = "A") & theme(legend.position = 'bottom')
dev.off()
