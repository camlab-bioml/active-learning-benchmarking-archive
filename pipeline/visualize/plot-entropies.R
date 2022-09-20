library(tidyverse)
library(RColorBrewer)
library(patchwork)

files <- snakemake@input[['f']]

## Create plots
create_boxplot <- function(vals, y_lab, colors){
  box <- vals %>% 
    mutate(no_cells_annotated = factor(no_cells_annotated)) %>% 
    ggplot(aes(x=no_cells_annotated, y=criterion_val, 
              group = no_cells_annotated, fill = no_cells_annotated)) + 
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = colors) +
    labs(y = y_lab) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
    
  box
}

if(snakemake@wildcards[['AL_type']] == 'Active-Learning_entropy'){
  ylab = "Entropy"
}else if(snakemake@wildcards[['AL_type']] == 'Active-Learning_maxp'){
  ylab = "Maxp"
}
plots <- lapply(files, function(x){
  vals <- read_tsv(x)
  col_numbers <- vals$no_cells_annotated %>% 
    unique %>% 
    length()
  colors <- colorRampPalette(brewer.pal(9, "Blues"))(col_numbers)

  create_boxplot(vals, ylab, colors)
})


pdf(snakemake@output[['box']], height = 9, width = 14)
  plots[[1]] / plots[[2]] / plots[[3]]
dev.off()
  