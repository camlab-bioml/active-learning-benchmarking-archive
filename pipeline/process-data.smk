
process_data_output = {
    'train_test_split': expand('data/{modality}/{modality}-{split}.rds', modality = modalities, split = data_splits),
    'scRNA_expression_df': expand('data/{modality}/{modality}-expression-df-{split}.tsv', modality = modalities, split = data_splits),
    'CyTOF_download': 'data/CyTOF/CyTOF-full.rds',
    'dataset_dimensionality': output + 'reports/dataset-dimensionality.html', 

    # Random subsets
    'sce_subset1': expand('data/{modality}/random/random-{modality}-set1.rds', modality = modalities),
    'sce_subset2': expand('data/{modality}/random/random-{modality}-set2.rds', modality = modalities),
    'sce_subset3': expand('data/{modality}/random/random-{modality}-set3.rds', modality = modalities),
    'gt_subset1': expand('data/scRNASeq/random/random-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-set1.tsv', modality = modalities),
    'gt_subset2': expand('data/scRNASeq/random/random-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-set2.tsv', modality = modalities),
    'gt_subset3': expand('data/scRNASeq/random/random-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-set3.tsv', modality = modalities),

    # Seurat clustering subsets
    'seu_sce': expand('data/{modality}/Seurat-clustering/Seurat-clustering-{modality}-max_dim-{pca_max_dim}-resolution-{res}.rds', 
        modality = ['scRNASeq'], pca_max_dim = [10, 12, 15], res = [0.3, 0.5, 0.7]),
    'gt_seu': expand('data/{modality}/Seurat-clustering/Seurat-clustering-{modality}-annotator-GroundTruth-max_dim-{pca_max_dim}-resolution-{res}-iterations_set-NA.tsv', 
        modality = ['scRNASeq'], pca_max_dim = [10, 12, 15], res = [0.3, 0.5, 0.7]),
}

rule split_datasets:
    input:
        rds = 'data/{modality}/{modality}-full.rds'
    resources:
        mem_mb=5000
    output:
        train = 'data/{modality}/{modality}-train.rds',
        test = 'data/{modality}/{modality}-test.rds',
        #random_train = 'data/{modality}/random/random-{modality}-500-cells.rds',
    script:
        'process-data/split-train-test.R'


rule process_data_for_random_forest:
    input:
        expression_sce = 'data/{modality}/{modality}-{split}.rds'
    output:
        expression_df = 'data/{modality}/{modality}-expression-df-{split}.tsv'
    script:
        'process-data/Create-expression-df.R'

rule download_CyTOF:
    output:
        Levine_CyTOF = 'data/CyTOF/CyTOF-full.rds'
    script:
        'process-data/download-and-process-CyTOF.R'


rule determine_dataset_dimensionality:
    input:
        markers = 'markers/{modality}.yml',
        training_rds = 'data/{modality}/{modality}-train.rds'
    output:
        html = output + 'reports/dataset-dimensionality.html'
    script:
        'process-data/determine-dataset-dimensionality.Rmd'

rule create_random_subsets:
    input:
        sce = 'data/{modality}/{modality}-train.rds'
    output:
        sce_subset1 = 'data/{modality}/random/random-{modality}-set1.rds',
        sce_subset2 = 'data/{modality}/random/random-{modality}-set2.rds',
        sce_subset3 = 'data/{modality}/random/random-{modality}-set3.rds',
        gt_subset1 = 'data/scRNASeq/random/random-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-set1.tsv',
        gt_subset2 = 'data/scRNASeq/random/random-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-set2.tsv',
        gt_subset3 = 'data/scRNASeq/random/random-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-set3.tsv'
    script:
        'process-data/select-random-subset.R'


rule create_clustering_subsets:
    input:
        seurat = output + 'cluster-and-interpret/{modality}/{modality}-Seurat-assignments-max_dim-{pca_max_dim}-resolution-{res}.tsv',
        sce = 'data/{modality}/{modality}-train.rds'
    output:
        ground_truth = 'data/{modality}/Seurat-clustering/Seurat-clustering-{modality}-annotator-GroundTruth-max_dim-{pca_max_dim}-resolution-{res}-iterations_set-NA.tsv',
        sce = 'data/{modality}/Seurat-clustering/Seurat-clustering-{modality}-max_dim-{pca_max_dim}-resolution-{res}.rds'
    script:
        'process-data/select-cluster-subset.R'


