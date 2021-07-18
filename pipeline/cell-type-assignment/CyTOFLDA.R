suppressPackageStartupMessages({
  library(scater)
  library(SingleCellExperiment)
  library(tidyverse)
  source(file.path("CyTOF-Linear-Classifier/", 'CyTOF_LDAtrain.R'))
  source(file.path("CyTOF-Linear-Classifier/", 'CyTOF_LDApredict.R'))
})

sce_train <- readRDS(snakemake@input[['training_rds']])

labs <- read_tsv(snakemake@input[['labels']])
labs <- labs %>% column_to_rownames("cell_id")

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
  prediction_params = "CyTOF-LDA",
  selection_procedure = snakemake@wildcards[['selection_procedure']],
  training_annotator = snakemake@wildcards[['annotator']],
  modality = 'CyTOF'
)

write_tsv(df_output, snakemake@output[['prediction']])