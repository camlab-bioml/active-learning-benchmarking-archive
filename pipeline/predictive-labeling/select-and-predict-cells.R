suppressPackageStartupMessages({
  library(caret)
  library(SingleCellExperiment)
  library(yaml)
  library(data.table)
  library(tibble)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggalluvial)
})
source("pipeline/whatsthatcell-helpers.R")

### [ READ IN DATA ] #####
markers <- read_yaml("markers/CyTOF.yml") %>% 
  unlist() %>% 
  unique()

sce <- readRDS("data/CyTOF/CyTOF-train.rds")
if(any(grepl("CD16_32", rownames(sce)))){
  rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
}

# Get selected set of cells to subset to for training
selected_cells <- read_tsv("data/CyTOF/Active-Learning/AL-batches-subset/Active-Learning-lowest_entropy-rand_sel-0-corr-0.1-CyTOF-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-500.tsv")


### [ PROCESS DATA ] #####
# create dataframe & add labels
df_expression <- load_scs(sce)
df_expression$cell_type <- sce$CellType

df_PCA <- select(df_expression, -X1) |> 
  as.matrix() |> 
  prcomp(center = TRUE, scale. = TRUE)

df_PCA <- df_PCA$x |> 
  as.data.frame()

df_PCA <- bind_cols(
  tibble(X1 = df_expression$X1),
  df_PCA[,1:20], 
  tibble(cell_type = sce$CellType)
)

# Separate into labelled and unlabeled 
labelled <- df_PCA %>% 
  filter(X1 %in% selected_cells$cell_id)

unlabeled <- df_PCA %>% 
  filter(!(X1 %in% selected_cells$cell_id))

### [ TRAIN AND PREDICT ] #####
# Train LR on labelled subset
multiNomModelFit <- train(cell_type ~ .,
                          data = select(labelled, -X1),
                          method = "multinom",
                          trace = FALSE)

# Predict probabilities
unlabeled_pred <- predict(multiNomModelFit,
                          select(unlabeled, -X1, -cell_type), 
                          type = "prob")

## Calculate entropy and max probability
entropy <- apply(unlabeled_pred, 1, calculate_entropy)

max_p_idx <- apply(unlabeled_pred, 1, which.max)
max_p <- lapply(seq_len(nrow(unlabeled_pred)), function(x){
  unlabeled_pred[x, max_p_idx[x]]
}) %>% unlist()

unlabeled_pred$entropy <- entropy
unlabeled_pred$max_p <- max_p

# Add predicted labels
unlabeled_pred$pred_cell_type <- predict(multiNomModelFit, unlabeled[, markers])

# Add cell ID
unlabeled_pred$cell_id <- unlabeled$X1

# Remove probabilities for each cell type
unlabeled_pred <- select(unlabeled_pred, cell_id, pred_cell_type, entropy, max_p)

# Calculate number of cells to select
additional_cell_num <- round(nrow(selected_cells) * 0.20)

entropy_cells <- unlabeled_pred %>% 
  arrange(entropy) %>% 
  slice_head(n = additional_cell_num) %>% 
  select(cell_id, pred_cell_type) %>% 
  mutate(selection_strat = "entropy") %>% 
  as_tibble()

maxp_cells <- unlabeled_pred %>% 
  arrange(-max_p) %>% 
  slice_head(n = additional_cell_num) %>% 
  select(cell_id, pred_cell_type) %>% 
  mutate(selection_strat = "maxp") %>% 
  as_tibble()

selected_cells <- bind_rows(entropy_cells, maxp_cells)

pdf(snakemake@output[[]])
  selected_cells %>% 
    left_join(select(unlabeled, X1, cell_type), by = c("cell_id" = "X1")) %>% 
    dplyr::rename("ground_truth" = "cell_type") %>% 
    pivot_longer(c(pred_cell_type, ground_truth), names_to = "type", values_to = "cell_type") %>% 
    ggplot(aes(x = type, stratum = cell_type, fill = cell_type, alluvium = cell_id)) +
    geom_stratum() +
    geom_flow() +
    scale_fill_manual(values = cell_type_colours()) +
    whatsthatcell_theme() +
    facet_wrap(~selection_strat) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

write_tsv(selected_cells, snakemake@output[[]])

