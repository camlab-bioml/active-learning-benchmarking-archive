
cell_type_predictions = {
    'random_forest_models': expand(output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl', 
        modality = modalities, selection_procedure = selection_procedures, annotator = annotators),
    'scmap_cluster': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-cluster-predictions.tsv',
        selection_procedure = selection_procedures, annotator = annotators),
    'scmap_sc': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-scmap-sc-predictions.tsv',
        selection_procedure = selection_procedures, annotator = annotators),
    'singleR': expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-annotation-{annotator}-singleR-predictions.tsv',
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
        'random-forest-train.py'

rule predict_random_forest:
    input:
        test = 'data/{modality}/{modality}-expression-df-test.tsv',
        model = output + 'models/random-forest-{modality}-trained-on-{selection_procedure}-by-{annotator}.pkl'
    output:
        predictions = output + 'rare-subtype-benchmarking/{modality}-{selection-procedure}-annotation-{annotator}-randomForest-predictions.tsv'
    script:
        'random-forest-test.py'