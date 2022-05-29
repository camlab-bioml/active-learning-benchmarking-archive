library(tidyverse)
library(yardstick)
source("pipeline/whatsthatcell-helpers.R")

### [ LOAD GROUND TRUTH ] #####
scRNAseq <- readRDS(snakemake@input[['sce']])
sce <- readRDS("data/CyTOF/CyTOF-test.rds")
expression_gt <- tibble(cell_id = colnames(sce),
                        annotated_cell_type = sce$CellType)

AL <- list.files("output/v3/rare-subtype-benchmarking/",
                 pattern = "CyTOF-Active-Learning", full.names = TRUE)

AL_files <- snakemake@input[['assignments']]
save.image('debug')
AL <- lapply(AL_files, read_tsv) %>% 
  bind_rows()

split_params <- str_split_fixed(AL$prediction_params, "-", 14)
AL$method <- paste0(split_params[,1], "-", split_params[,2])
AL$cell_num <- split_params[,10]
AL$randomSelection <- split_params[,12]
AL$corrupted <- split_params[,14]

AL <- left_join(AL, expression_gt)


acc <- AL %>% 
  filter(predicted_cell_type != "unassigned") %>% 
  group_by(method, selection_procedure, cell_num, randomSelection, corrupted) %>% 
  acc_wrap()

pdf(snakemake@output[['pdf']], width = 12, height = 4)
  acc %>% 
    mutate(cell_num = factor(as.integer(cell_num))) %>% 
    filter(method == "CyTOF-LDA") %>% 
    ggplot(aes(x = cell_num, y = .estimate, color = selection_procedure, group = selection_procedure)) +
    geom_point() +
    geom_line(aes(group = selection_procedure)) +
    scale_x_discrete(breaks = acc$cell_num) +
    ylim(0,1) +
    theme_bw() +
    facet_grid(randomSelection + corrupted~.metric)
dev.off()
