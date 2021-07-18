library(HDCytoData)
library(SingleCellExperiment)

# Download dataset
CyTOF <- Levine_32dim_SE()

# remove non-celltype markers & unassigned cells
assigned <- CyTOF[rowData(CyTOF)$population_id != "unassigned", 
                  colData(CyTOF)$marker_class == "type"]

# Convert to single cell experiment & normalize
expression <- assay(assigned)
cell_types <- as.character(rowData(assigned)$population_id)

# Normalize
expression <- asinh(expression / 5)

sce <- SingleCellExperiment(list(logcounts = t(expression)),
                            colData = DataFrame(cell_type = cell_types))
sce$cell_type <- cell_types
colnames(sce) <- paste0("Levine_32_", seq(1:ncol(sce)))

saveRDS(sce, snakemake@output[['Levine_CyTOF']])



### QUESTIONS
# - should I be normalizing these data?
# - Can I normalize after removal of unassigned?