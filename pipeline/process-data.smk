

process_data_output = {
    'train_test_split': expand('data/{modality}/{modality}-{split}-seed-{s}.rds', modality = modalities, split = data_splits, s = train_test_seeds),
    'scRNA_expression_df': expand('data/{modality}/{modality}-expression-df-{split}-seed-{s}.tsv', modality = modalities, split = data_splits, s = train_test_seeds),
    'CyTOF_download': 'data/CyTOF/CyTOF-full.rds',
    'dataset_dimensionality': expand(output + 'reports/{modality}-dataset-dimensionality-seed-{s}.html', modality = modalities, s = train_test_seeds),

    # Random subsets
    'sce_subset1': expand(output + 'data/{modality}/random/random-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-{cell_num}_cells-seed-{s}.rds', modality = modalities, corrupt = [0], cell_num = cell_numbers, set = random_sets, s = train_test_seeds),
    'gt_subset1': expand(output + 'data/{modality}/random/random-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations-{cell_num}_cells-seed-{s}.tsv', modality = modalities, corrupt = [0], cell_num = cell_numbers, s = train_test_seeds),
    
    # Seurat clustering subsets
    'seu_sce': expand(output + 'data/{modality}/Seurat-clustering/Seurat-clustering-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-knn_neighbors-{neighbors}-resolution-{res}-{cell_num}_cells-seed-{s}.rds', 
        modality = modalities, neighbors = Seurat_neighbors, res = Seurat_resolution, corrupt = [0], cell_num = cell_numbers, s = train_test_seeds),
    'gt_seu': expand(output + 'data/{modality}/Seurat-clustering/Seurat-clustering-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-NA-{cell_num}_cells-seed-{s}.tsv', 
        modality = modalities, neighbors = Seurat_neighbors, res = Seurat_resolution, corrupt = [0], cell_num = cell_numbers, s = train_test_seeds),
}

rule split_datasets:
    input:
        rds = 'data/{modality}/{modality}-full.rds'
    resources:
        mem_mb=5000
    output:
        train = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test = 'data/{modality}/{modality}-test-seed-{s}.rds',
    script:
        'process-data/split-train-test.R'

rule process_data_for_random_forest:
    input:
        expression_sce = 'data/{modality}/{modality}-{split}-seed-{s}.rds'
    resources:
        mem_mb=5000
    output:
        expression_df = 'data/{modality}/{modality}-expression-df-{split}-seed-{s}.tsv'
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
        training_rds = 'data/{modality}/{modality}-train-seed-{s}.rds'
    resources:
        mem_mb=5000
    output:
        html = output + 'reports/{modality}-dataset-dimensionality-seed-{s}.html'
    script:
        'process-data/determine-dataset-dimensionality.Rmd'

rule create_random_subsets:
    input:
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    output:
        sce_subset1 = output + 'data/{modality}/random/random-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-{cell_num}_cells-seed-{s}.rds',
        gt_subset1 = output + 'data/{modality}/random/random-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations-{cell_num}_cells-seed-{s}.tsv',
    script:
        'process-data/select-random-subset.R'

rule create_clustering_subsets:
    input:
        seurat = output + 'cluster-and-interpret/{modality}/{modality}-Seurat-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    output:
        ground_truth = output + 'data/{modality}/Seurat-clustering/Seurat-clustering-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-NA-{cell_num}_cells-seed-{s}.tsv',
        sce = output + 'data/{modality}/Seurat-clustering/Seurat-clustering-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{modality}-knn_neighbors-{neighbors}-resolution-{res}-{cell_num}_cells-seed-{s}.rds'
    script:
        'process-data/select-cluster-subset.R'
