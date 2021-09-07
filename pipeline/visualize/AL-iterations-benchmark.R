library(tidyverse)
library(yardstick)
source("pipeline/whatsthatcell-helpers.R")

### [ LOAD GROUND TRUTH ] #####
scRNAseq <- readRDS(snakemake@input[['sce']])
scRNAseq_gt <- tibble(cell_id = colnames(scRNAseq),
                      annotated_cell_type = scRNAseq$CellType)

# AL <- list.files("output/v1/rare-subtype-benchmarking/", 
#                  pattern = "scRNASeq-Active-Learning", full.names = TRUE)

AL_files <- snakemake@input[['assignments']]
save.image('debug')
AL <- lapply(AL_files, read_tsv) %>% 
  bind_rows()

split_params <- str_split_fixed(AL$prediction_params, "-", 4)
AL$method <- paste0(split_params[,1], "-", split_params[,2])
AL$iteration <- split_params[,4]

AL <- left_join(AL, scRNAseq_gt)


acc <- AL %>% 
  filter(predicted_cell_type != "unassigned") %>% 
  group_by(method, iteration) %>% 
  acc_wrap()

pdf(snakemake@output[['pdf']], width = 12, height = 4)
  acc %>% 
    mutate(iteration = factor(as.integer(iteration))) %>% 
    ggplot(aes(x = iteration, y = .estimate, group = iteration)) +
    geom_point() +
    scale_x_discrete(breaks = acc$iteration) +
    theme_bw() +
    facet_grid(method~.metric)
dev.off()
