
process_data_output = {
    'train_test_split': expand('data/{modality}/{modality}-{split}.rds', modality = modalities, split = data_splits),
    'random_subset': expand('data/{modality}/random/random-{modality}-500-cells.rds', modality = modalities),
    'scRNA_expression_df': expand('data/{modality}/{modality}-expression-df-{split}.tsv', modality = modalities, split = data_splits)
}

rule split_datasets:
    input:
        rds = 'data/scRNASeq/scRNASeq-full.rds'
    resources:
        mem_mb=5000
    output:
        train = 'data/scRNASeq/scRNASeq-train.rds',
        test = 'data/scRNASeq/scRNASeq-test.rds',
        random_train = 'data/scRNASeq/random/random-scRNASeq-500-cells.rds',
    script:
        'process-data/split-train-test.R'


rule process_data_for_random_forest:
    input:
        expression_sce = 'data/{modality}/{modality}-full.rds'
    output:
        expression_df = 'data/{modality}/{modality}-expression-df-{split}.tsv'
    script:
        'process-data/Create-expression-df.R'
