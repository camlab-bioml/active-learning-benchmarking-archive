library(tidyverse)
library(RColorBrewer)

entropy <- read_tsv(snakemake@input[['entropy']])

col_numbers <- entropy$no_cells_annotated %>% 
  unique %>% 
  length
colors <- colorRampPalette(brewer.pal(9, "Blues"))(col_numbers) %>% 
  rev()

pdf(snakemake@output[['box']], height = 4, width = 8)
entropy %>% 
  mutate(no_cells_annotated = factor(no_cells_annotated)) %>% 
  ggplot(aes(x=no_cells_annotated, y=entropy, 
                    group = no_cells_annotated, fill = no_cells_annotated)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x = "Number of cells annotated", y = "Entropy") +
  theme(legend.position = "none")
dev.off()
  

pdf(snakemake@output[['hist']], height = 20, width = 10)
ggplot(entropy, aes(x = entropy)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  facet_wrap(~no_cells_annotated, ncol = 10)
dev.off()