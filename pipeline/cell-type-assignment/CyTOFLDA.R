suppressPackageStartupMessages({
  library(scater)
  library(SingleCellExperiment)
  library(tidyverse)
  source(file.path("CyTOF-Linear-Classifier/", 'CyTOF_LDAtrain.R'))
  source(file.path("CyTOF-Linear-Classifier/", 'CyTOF_LDApredict.R'))
})
library(devtools)
devtools::load_all("../taproom/")

sce_train <- readRDS("data/training/zurich_subset1k.rds")
labs <- read_csv("data/zurich1_astir_assignments.csv")
cell_id <- labs$X1
labs$X1 <- NULL
labs$cell_type <- get_celltypes(labs)
labs <- as.data.frame(labs)
rownames(labs) <- cell_id


labs <- select(labs, cell_type)

sce_train$cell_type <- labs[colnames(sce_train), 'cell_type']

sce_annotate <- sce_train[,1:500]
sce_train <- sce_train[,501:1000]


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

predictions <- unlist(Predictions)

df_output <- tibble(
  cell_id = data_annotate$cell_id,
  cell_type = Predictions,
  annotator = args$annotator,
  cohort = args$cohort,
  method = args$method
)

write_tsv(df_output, args$output_assignments)