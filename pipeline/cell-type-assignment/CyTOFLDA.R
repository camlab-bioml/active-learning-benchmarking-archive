suppressPackageStartupMessages({
  library(scater)
  library(SingleCellExperiment)
  library(tidyverse)
  library(caret)
  source(file.path("utils/", 'CyTOF_LDAtrain.R'))
  source(file.path("utils/", 'CyTOF_LDApredict.R'))
})
sce_train <- readRDS(snakemake@input[['training_rds']])

labs <- read_tsv(snakemake@input[['labels']])
labs <- labs %>% column_to_rownames("cell_id")

# # For active learning, make sure the cells are in the order they were selected in
# if(snakemake@wildcards[['selection_procedure']] == "Active-Learning"){
#   labs <- labs %>% arrange(iteration)
# }
# print('this owrks')
# if(snakemake@wildcards[['selection_procedure']] == "Seurat-clustering"){
#   # For seurat clustering - make sure I still have an even number of cell types
#   if(snakemake@wildcards[['cell_num']] == "500"){
#     prop <- 1
#     print('mayby')
#   }else{
#     print('mabye2')
#     prop <- as.integer(snakemake@wildcards[['cell_num']]) / nrow(labs)
#   }
#   print('maybe 3')
#   trainIndex <- createDataPartition(labs$cell_type, p = prop, 
#                                     list = FALSE, 
#                                     times = 1)
#   labs <- labs[trainIndex, ]
# }else{
#   print('maybe 5')
#   # Subset to the number of cells that are of interest
#   labs <- labs[1:as.integer(snakemake@wildcards[['cell_num']]),]
# }

sce_train <- sce_train[,colnames(sce_train) %in% rownames(labs)]
sce_train$cell_type <- labs[colnames(sce_train), 'cell_type']

sce_annotate <- read_rds(snakemake@input[['annotation_rds']])

data_train <- as.data.frame(t(logcounts(sce_train)))
data_train$cell_type <- sce_train$cell_type

data_annotate <- as.data.frame(t(logcounts(sce_annotate)))

train_dir <- file.path(tempdir(), "train")
annotate_dir <- file.path(tempdir(), "annotate")

dir.create(train_dir)
dir.create(annotate_dir)

write.table(data_train, file = file.path(train_dir, "train.csv"), 
            col.names = FALSE, row.names = FALSE, sep = ',')
write.table(data_annotate, file = file.path(annotate_dir, "annotate.csv"),
            col.names = FALSE, row.names = FALSE, sep = ',')

LDA.Model <- CyTOF_LDAtrain(TrainingSamplesExt = train_dir, TrainingLabelsExt = '',
                            mode = 'CSV', RelevantMarkers = seq_len(nrow(sce_train)),
                            LabelIndex = ncol(data_train), Transformation = FALSE)

predictions <- CyTOF_LDApredict(LDA.Model, TestingSamplesExt = annotate_dir,
                                mode = 'CSV', RejectionThreshold = 0)

predictions <- unlist(predictions)

df_output <- tibble(
  cell_id = rownames(data_annotate),
  predicted_cell_type = predictions,
  prediction_params = paste0("CyTOF-LDA-iterations_set-", snakemake@wildcards[['set']], "-knn-", 
                             snakemake@wildcards[['neighbors']], "-res-", snakemake@wildcards[['res']], '-cell_numbers-', snakemake@wildcards[['cell_num']],
                             '-randomSelection-', snakemake@wildcards[['rand']], '-corrupted-', snakemake@wildcards[['corrupt']]),
  selection_procedure = paste0(snakemake@wildcards[['selection_procedure']], '-strategy-', snakemake@wildcards[['strat']]),
  training_annotator = snakemake@wildcards[['annotator']],
  modality = 'CyTOF'
)

write_tsv(df_output, snakemake@output[['prediction']])