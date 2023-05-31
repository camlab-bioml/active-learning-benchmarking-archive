suppressPackageStartupMessages({
  library(scater)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(magick)
})
source("pipeline/whatsthatcell-helpers.R")


### [ FIGURE 1A - BENCHMARKING OVERVIEW ] #####
schematic <- image_read("illustrator-figures/benchmarking-schematic.ai") |> 
  image_ggplot()

### [ FIGURE 1B - DATASET COMPOSOTION ] #####
set.seed(42)

CyTOF <- readRDS("data/CyTOF/CyTOF-full.rds")
CyTOF <- runTSNE(CyTOF)
scRNA <- readRDS("data/scRNASeq/scRNASeq-full.rds")
scRNA <- runTSNE(scRNA)
snRNA <- readRDS("data/snRNASeq/snRNASeq-full.rds")
snRNA <- runTSNE(snRNA)


plot_dim_red <- function(sce, mod, include_axis = FALSE,
                         a1_start = 0, a1_end = 0, a1_y = 0,
                         a2_start = 0, a2_end = 0, a2_x = 0, l_nrow = 4){
  if(!("CellType" %in% names(colData(sce)))){
    sce$CellType <- sce$cell_type
  }
  p <- plotTSNE(sce, colour_by = 'CellType') +
    scale_color_manual(values = cell_type_colours(mod, FALSE)) +
    coord_fixed() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank()) + 
    guides(color = guide_legend(nrow = l_nrow, byrow = TRUE))
  
  if(include_axis){
    p <- p + geom_segment(x = a1_start, xend = a1_end, y = a1_y, yend = a1_y,
                     arrow = arrow(length = unit(0.2, "cm")),
                     lineend = "round", linejoin = "mitre", size = 1) +
      geom_segment(x = a2_x, xend = a2_x, y = a2_start, yend = a2_end,
                   arrow = arrow(length = unit(0.2, "cm")),
                   lineend = "round", linejoin = "mitre", size = 1) +
      theme(axis.title.x = element_text(hjust = 0.035),
            axis.title.y = element_text(hjust = 0.035, vjust = -1.5))
  }else{
    p <- p +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
  }
  
  p
}



### CYTOF
cytof_tsne <- plot_dim_red(CyTOF, "CyTOF")
cytof_bar <- CyTOF$cell_type |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 400),
            position = position_stack(vjust = 1.03)) +
  scale_fill_manual(values = cell_type_colours("CyTOF", FALSE)) +
  ylim(0, 2800) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

cytof_plot <- (wrap_elements(full = cytof_tsne, ignore_tag = TRUE) & labs(title = "CyTOF")) /
  cytof_bar + plot_layout(heights = c(3,1))



### scRNASeq
scrna_tsne <- plot_dim_red(scRNA, "scRNASeq", TRUE,
                           a1_start = -41, a1_end = -33, a1_y = -47,
                           a2_start = -47, a2_end = -38, a2_x = -41, 
                           l_nrow = 5)

scrna_bar <- scRNA$CellType |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 50),
            position = position_stack(vjust = 1.01)) +
  scale_fill_manual(values = cell_type_colours("scRNASeq", FALSE)) +
  ylim(0, 1650) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


scrna_plot <- ((wrap_elements(full = scrna_tsne, ignore_tag = TRUE) & labs(title = "scRNASeq")) /
  scrna_bar) + 
  plot_layout(heights = c(3,1))


### snRNASeq
snrna_tsne <- plot_dim_red(snRNA, "snRNASeq")

snrna_bar <- snRNA$cell_type |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 400),
            position = position_stack(vjust = 1.03)) +
  scale_fill_manual(values = cell_type_colours("snRNASeq", FALSE)) +
  ylim(0, 3500) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        #legend.position = "none",
        axis.ticks.x = element_blank())


snrna_plot <- ((wrap_elements(full = snrna_tsne, ignore_tag = TRUE) & labs(title = "snRNASeq")) / 
                 snrna_bar) + 
  plot_layout(heights = c(3,1))

pdf("output/v8/paper-figures/figure1.pdf", height = 12.5, width = 17.5)
  wrap_elements((scrna_plot | snrna_plot | cytof_plot) + plot_layout(widths = c(1, 1.2, 1))) /
    wrap_elements(schematic) +
    plot_layout(heights = c(2, 0.7)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 22))
dev.off()


