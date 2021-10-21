
selection_expansion_dict = {
    'Seurat-clustering': {'neighbors': Seurat_neighbors,
                        'res': Seurat_resolution,
                        'set': ['NA']},
    'random': {'neighbors': ['NA'],
               'res': ['NA'],
               'set': ['set1', 'set2', 'set3']},
    'Active-Learning': {'neighbors': ['NA'],
               'res': ['NA'],
               'set': list(range(15))}
}

scRNASeq_methods = ['scmap-cluster', 'scmap-sc', 'singleR']

scRNA_predictions = []
scRNA = [expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-{method}-predictions.tsv',
                        selection_procedure = [select], annotator = ['GroundTruth'], 
                        neighbors = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        set = selection_expansion_dict[select]['set'],
                        method = scRNASeq_methods) 
                        for select in selection_expansion_dict.keys()]
for element in scRNA:
    scRNA_predictions.extend(element)


cell_type_predictions = {
    'scRNASeq': scRNA_predictions
    #'random_forest_models': expand(output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl', 
    #    modality = modalities, selection_procedure = selection_procedures, annotator = annotators),
    #'random_forest_prediction': expand(output + 'rare-subtype-benchmarking/{modality}-{selection_procedure}-annotation-{annotator}-randomForest-predictions.tsv',
    #    modality = modalities, selection_procedure = selection_procedures, annotator = annotators),
    # 'scmap_cluster': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-cluster-predictions.tsv',
    #     selection_procedure = selection_procedures, annotator = annotators),
    # 'scmap_cluster_AL': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-cluster-predictions.tsv',
    #     selection_procedure = ['Active-Learning'], annotator = AL_annotator),
    # 'scmap_sc': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-sc-predictions.tsv',
    #     selection_procedure = selection_procedures, annotator = annotators),
    # 'scmap_sc': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-sc-predictions.tsv',
    #     selection_procedure = ['Active-Learning'], annotator = AL_annotator),
    # 'singleR': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-singleR-predictions.tsv',
    #     selection_procedure = selection_procedures, annotator = annotators),
    # 'singleR': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-singleR-predictions.tsv',
    #     selection_procedure = ['Active-Learning'], annotator = AL_annotator),
    # 'Seurat': expand(output + 'cluster-and-interpret/{modality}/{modality}-Seurat-assignments-max_dim-{max_pca_dim}-resolution-{res}.tsv', 
    #     modality = ['scRNASeq'], max_pca_dim = [10, 12, 15], res = [0.3, 0.5, 0.7]),
    # 'CyTOF_LDA': expand(output + 'rare-subtype-benchmarking/CyTOF-{selection_procedure}-annotation-{annotator}-CyTOF-LDA-predictions.tsv', 
    #     selection_procedure = selection_procedures, annotator = annotators)
}

rule train_and_predict_scmap:
    input:
        annotation = 'data/scRNASeq/{selection_procedure}/{selection_procedure}-scRNASeq-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}.tsv',
        train_data = 'data/scRNASeq/scRNASeq-train.rds',
        test_data = 'data/scRNASeq/scRNASeq-test.rds'
    output:
        cluster_predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-scmap-cluster-predictions.tsv',
        sc_predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-scmap-sc-predictions.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule train_and_predict_singleR:
    input:
        annotation = 'data/scRNASeq/{selection_procedure}/{selection_procedure}-scRNASeq-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}.tsv',
        train_data = 'data/scRNASeq/scRNASeq-train.rds',
        test_data = 'data/scRNASeq/scRNASeq-test.rds'
    output:
        predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-singleR-predictions.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule train_random_forest:
    input:
        train = 'data/{modality}/{modality}-expression-df-train.tsv',
        annotation = 'data/{modality}/{selection_procedure}/{selection_procedure}-annotation-{annotator}.tsv'
    output:
        model = output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl'
    resources:
        mem_mb=20000
    log:
        output + 'logs/cell-type-predictions/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule predict_random_forest:
    input:
        test = 'data/{modality}/{modality}-expression-df-test.tsv',
        model = output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl'
    resources:
        mem_mb=5000
    output:
        predictions = output + 'rare-subtype-benchmarking/{modality}-{selection_procedure}-annotation-{annotator}-randomForest-predictions.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule Seurat_clustering:
    input:
        training_rds = 'data/{modality}/{modality}-train.rds',
        markers = 'markers/{modality}.yml'
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/{modality}-Seurat-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/{modality}-Seurat-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}.pdf',
    output:
        cluster_umap_pdf = output + 'figures/{modality}-Seurat-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        cell_type_umap_pdf = output + 'figures/{modality}-Seurat-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        assignments = output + 'cluster-and-interpret/{modality}/{modality}-Seurat-assignments-knn_neighbors-{neighbors}-resolution-{res}.tsv',
        diagnostics = expand(output + 'figures/diagnostics/{{modality}}-Seurat-{cell_types}-{pn}-knn_neighbors-{{neighbors}}-resolution-{{res}}.pdf', 
            cell_types = scRNA_cell_types_clean, pn = ['positive']),
        ground_truth_umap_pdf = output + 'figures/{modality}-Seurat-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf'
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
        
