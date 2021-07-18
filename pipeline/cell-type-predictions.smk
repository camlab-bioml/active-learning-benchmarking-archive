
cell_type_predictions = {
    'random_forest_models': expand(output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl', 
        modality = modalities, selection_procedure = selection_procedures, annotator = annotators),
    'random_forest_prediction': expand(output + 'rare-subtype-benchmarking/{modality}-{selection_procedure}-annotation-{annotator}-randomForest-predictions.tsv',
        modality = modalities, selection_procedure = selection_procedures, annotator = annotators),
    'scmap_cluster': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-cluster-predictions.tsv',
        selection_procedure = selection_procedures, annotator = annotators),
    'scmap_sc': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-sc-predictions.tsv',
        selection_procedure = selection_procedures, annotator = annotators),
    'singleR': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-singleR-predictions.tsv',
        selection_procedure = selection_procedures, annotator = annotators),
    'Seurat': expand(output + 'cluster-and-interpret/Seurat-assignments-max_dim-{max_pca_dim}-resolution-{res}.tsv', max_pca_dim = [15], res = [0.5]),
    'CyTOF_LDA': expand(output + 'rare-subtype-benchmarking/CyTOF-{selection_procedure}-annotation-{annotator}-CyTOF-LDA-predictions.tsv', 
        selection_procedure = selection_procedures, annotator = annotators)
}

rule train_and_predict_scmap:
    input:
        annotation = 'data/scRNASeq/{selection_procedure}/{selection_procedure}-annotation-{annotator}.tsv',
        train_data = 'data/scRNASeq/scRNASeq-train.rds',
        test_data = 'data/scRNASeq/scRNASeq-test.rds'
    output:
        cluster_predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-cluster-predictions.tsv',
        sc_predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-sc-predictions.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule train_and_predict_singleR:
    input:
        annotation = 'data/scRNASeq/{selection_procedure}/{selection_procedure}-annotation-{annotator}.tsv',
        train_data = 'data/scRNASeq/scRNASeq-train.rds',
        test_data = 'data/scRNASeq/scRNASeq-test.rds'
    output:
        predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-singleR-predictions.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule train_random_forest:
    input:
        train = 'data/{modality}/{modality}-expression-df-train.tsv',
        annotation = 'data/{modality}/{selection_procedure}/{selection_procedure}-annotation-{annotator}.tsv'
    output:
        model = output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl'
    log:
        output + 'logs/cell-type-predictions/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule predict_random_forest:
    input:
        test = 'data/{modality}/{modality}-expression-df-test.tsv',
        model = output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl'
    output:
        predictions = output + 'rare-subtype-benchmarking/{modality}-{selection_procedure}-annotation-{annotator}-randomForest-predictions.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule Seurat_clustering:
    input:
        training_rds = 'data/scRNASeq/scRNASeq-train.rds',
        markers = 'markers/scRNA.yml'
    params:
        max_dim = 15,
        resolution = 0.5,
        positive_markers_diagnostic = output + 'figures/diagnostics/Seurat-[cell_types]-positive-max_dim-{max_pca_dim}-resolution-{res}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/Seurat-[cell_types]-negative-max_dim-{max_pca_dim}-resolution-{res}.pdf',
    output:
        cluster_umap_pdf = output + 'figures/Seurat-cluster-assignment-umap-max_dim-{max_pca_dim}-resolution-{res}.pdf',
        cell_type_umap_pdf = output + 'figures/Seurat-cell-assignment-umap-max_dim-{max_pca_dim}-resolution-{res}.pdf',
        assignments = output + 'cluster-and-interpret/Seurat-assignments-max_dim-{max_pca_dim}-resolution-{res}.tsv',
        diagnostics = expand(output + 'figures/diagnostics/Seurat-{cell_types}-{pn}-max_dim-{{max_pca_dim}}-resolution-{{res}}.pdf', 
            cell_types = scRNA_cell_types_clean, pn = ['positive']),
        ground_truth_umap_pdf = output + 'figures/Seurat-ground-truth-umap-max_dim-{max_pca_dim}-resolution-{res}.pdf'
    script:
        'cell-type-assignment/Seurat.R'


rule CyTOF_LDA:
    input:
        training_rds = 'data/CyTOF/CyTOF-train.rds',
        annotation_rds = 'data/CyTOF/CyTOF-test.rds',
        labels = 'data/CyTOF/{selection_procedure}/{selection_procedure}-annotation-{annotator}.tsv'
    output:
        prediction = output + 'rare-subtype-benchmarking/CyTOF-{selection_procedure}-annotation-{annotator}-CyTOF-LDA-predictions.tsv'
    script:
        'cell-type-assignment/CyTOFLDA.R'
        
