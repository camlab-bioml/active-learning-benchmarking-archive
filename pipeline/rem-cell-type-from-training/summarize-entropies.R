suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})
source("pipeline/whatsthatcell-helpers.R")

entropies <- lapply(snakemake@input$entropy_files, function(x){
  read_tsv(x) |> 
    filter(no_cells_annotated == 20) |> 
    mutate(rem_cell_type = gsub('AL_alg-(.*?)-strat-(.*?)-init-(.*?)-rem_celltype-',
                                "", params),
           rem_cell_type = gsub("-seed-[0-9]", "", rem_cell_type),
           AL_alg = gsub("AL_alg-", "", params),
           AL_alg = gsub("-strat-.*", "", AL_alg),
           strat = gsub("AL_alg-(.*?)-strat-", "", params),
           strat = gsub("-init-.*", "", strat),
           initial = gsub("(.*?)-init-", "", params),
           initial = gsub("-rem_celltype.*", "", initial),
           seed = gsub(".*-seed-", "", params)) |> 
    select(-c(criterion, no_cells_annotated, params, comp))
}) |> bind_rows()

pdf(snakemake@output$entropies, height = 10, width = 8)
  mutate(entropies, gt_rem = rem_cell_type == gt_cell_type) |> 
    ggplot(aes(x = as.factor(num_missing_cells), y = criterion_val, fill = gt_rem)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#ffa600", "#9763ff")) +
    labs(x = "Number of missing cells", y = "Scaled entropy", 
         fill = "Cell type\nremoved?", title = snakemake@wildcards$modality) +
    facet_grid(strat ~ AL_alg + initial) +
    whatsthatcell_theme()
dev.off()

