
selection_expansion_dict = {
    'Seurat-clustering': {'neighbors': Seurat_neighbors,
                        'res': Seurat_resolution,
                        'set': ['NA'],
                        'strategy': 'NA',
                        'random_selection': 'NA',
                        'corruption': [0]},
    'random': {'neighbors': ['NA'],
               'res': ['NA'],
               'set': random_sets,
               'strategy': 'NA',
               'random_selection': 'NA',
               'corruption': [0]},
    'Active-Learning_entropy': {'neighbors': ['NA'],
               'res': ['NA'],
               'set': ['NA'],
               'strategy': ['0.25_quant_entropy', '0.5_quant_entropy', '0.75_quant_entropy', 'highest_entropy'],
               'random_selection': [0, 0.25, 0.5, 0.75],
               'corruption': corruption_percentages},
    'Active-Learning_maxp': {'neighbors': ['NA'],
               'res': ['NA'],
               'set': ['NA'],
               'strategy': ['0.25_quant_maxp', '0.5_quant_maxp', '0.75_quant_maxp', 'lowest_maxp'],
               'random_selection': [0, 0.25, 0.5, 0.75],
               'corruption': corruption_percentages}
}

scRNA_predictions = []
scRNA = [expand(output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-{method}-predictions-{cell_num}-cells.tsv',
                        selection_procedure = [select], 
                        annotator = ['GroundTruth'], 
                        strat = selection_expansion_dict[select]['strategy'],
                        rand = selection_expansion_dict[select]['random_selection'],
                        corrupt = selection_expansion_dict[select]['corruption'],
                        neighbors = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        set = selection_expansion_dict[select]['set'],
                        method = evaluation_methods_dict['scRNASeq'],
                        cell_num = cell_numbers) 
                        for select in selection_expansion_dict.keys()]
for element in scRNA:
    scRNA_predictions.extend(element)

CyTOF_predictions = []
CyTOF = [expand(output + 'rare-subtype-benchmarking/CyTOF-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-{method}-predictions-{cell_num}-cells.tsv',
                        selection_procedure = [select], 
                        annotator = ['GroundTruth'], 
                        strat = selection_expansion_dict[select]['strategy'],
                        rand = selection_expansion_dict[select]['random_selection'],
                        corrupt = selection_expansion_dict[select]['corruption'],
                        neighbors = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        set = selection_expansion_dict[select]['set'],
                        method = evaluation_methods_dict['CyTOF'],
                        cell_num = cell_numbers) 
                        for select in selection_expansion_dict.keys()]
for element in CyTOF:
    CyTOF_predictions.extend(element)

cell_type_predictions = {
    'scRNASeq': scRNA_predictions,
    'CyTOF': CyTOF_predictions
}

def get_labels(procedure, mod):
    if procedure == "Active-Learning_entropy" or procedure == "Active-Learning_maxp":
        path = 'data/{modality}/{{selection_procedure}}/AL-batches-subset/{{selection_procedure}}-{{strat}}-rand_sel-{{rand}}-corr-{{corrupt}}-{modality}-annotator-{{annotator}}-knn_neighbors-{{neighbors}}-resolution-{{res}}-iterations_set-NA-{{cell_num}}_cells.tsv'
    else:
        path = 'data/{modality}/{{selection_procedure}}/{{selection_procedure}}-{{strat}}-rand_sel-{{rand}}-corr-{{corrupt}}-{modality}-annotator-{{annotator}}-knn_neighbors-{{neighbors}}-resolution-{{res}}-iterations_set-{{set}}-{{cell_num}}_cells.tsv'
    
    path = expand(path, modality = mod)
    return path

rule train_and_predict_scmap:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, 'scRNASeq'),
        train_data = 'data/scRNASeq/scRNASeq-train.rds',
        test_data = 'data/scRNASeq/scRNASeq-test.rds'
    output:
        cluster_predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-scmap-cluster-predictions-{cell_num}-cells.tsv',
        sc_predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-scmap-sc-predictions-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule train_and_predict_singleR:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, 'scRNASeq'),
        train_data = 'data/scRNASeq/scRNASeq-train.rds',
        test_data = 'data/scRNASeq/scRNASeq-test.rds'
    output:
        predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-singleR-predictions-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule train_random_forest_cytof:
    input:
        train = 'data/CyTOF/CyTOF-expression-df-train.tsv',
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, 'CyTOF'),
    params:
        modality = 'CyTOF'
    output:
        model = output + 'models/random-forest-CyTOF-trained-on-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-by-{annotator}-{cell_num}-cells.pkl'
    resources:
        mem_mb=20000
    log:
        output + 'logs/cell-type-predictions/random-forest-CyTOF-trained-on-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-by-{annotator}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule train_random_forest_scRNA:
    input:
        train = 'data/scRNASeq/scRNASeq-expression-df-train.tsv',
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, 'scRNASeq'),
    params:
        modality = 'scRNASeq'
    output:
        model = output + 'models/random-forest-scRNASeq-trained-on-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-by-{annotator}-{cell_num}-cells.pkl'
    resources:
        mem_mb=20000
    log:
        output + 'logs/cell-type-predictions/random-forest-scRNASeq-trained-on-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-by-{annotator}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/random-forest-train.py'


rule predict_random_forest_cytof:
    input:
        test = 'data/CyTOF/CyTOF-expression-df-test.tsv',
        model = output + 'models/random-forest-CyTOF-trained-on-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-by-{annotator}-{cell_num}-cells.pkl'
    params:
        modality = 'CyTOF'
    resources:
        mem_mb=5000
    output:
        predictions = output + 'rare-subtype-benchmarking/CyTOF-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-Random-Forest-predictions-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule predict_random_forest_scrna:
    input:
        test = 'data/scRNASeq/scRNASeq-expression-df-test.tsv',
        model = output + 'models/random-forest-scRNASeq-trained-on-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-by-{annotator}-{cell_num}-cells.pkl'
    params:
        modality = 'scRNASeq'
    resources:
        mem_mb=5000
    output:
        predictions = output + 'rare-subtype-benchmarking/scRNASeq-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-Random-Forest-predictions-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule Seurat_clustering_RNASeq:
    input:
        training_rds = 'data/scRNASeq/scRNASeq-train.rds',
        markers = 'markers/scRNASeq.yml'
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/scRNASeq-Seurat-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/scRNASeq-Seurat-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        mod = "scRNASeq"
    output:
        cluster_umap_pdf = output + 'figures/scRNASeq-Seurat-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        cell_type_umap_pdf = output + 'figures/scRNASeq-Seurat-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        assignments = output + 'cluster-and-interpret/scRNASeq/scRNASeq-Seurat-assignments-knn_neighbors-{neighbors}-resolution-{res}.tsv',
        diagnostics = expand(output + 'figures/diagnostics/scRNASeq-Seurat-{cell_types}-{pn}-knn_neighbors-{{neighbors}}-resolution-{{res}}.pdf', 
                             cell_types = all_cell_types['scRNASeq'], pn = ['positive']),
        ground_truth_umap_pdf = output + 'figures/scRNASeq-Seurat-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_CyTOF:
    input:
        training_rds = 'data/CyTOF/CyTOF-train.rds',
        markers = 'markers/CyTOF.yml'
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/CyTOF-Seurat-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/CyTOF-Seurat-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        mod = "CyTOF"
    output:
        cluster_umap_pdf = output + 'figures/CyTOF-Seurat-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        cell_type_umap_pdf = output + 'figures/CyTOF-Seurat-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf',
        assignments = output + 'cluster-and-interpret/CyTOF/CyTOF-Seurat-assignments-knn_neighbors-{neighbors}-resolution-{res}.tsv',
        diagnostics = expand(output + 'figures/diagnostics/CyTOF-Seurat-{cell_types}-{pn}-knn_neighbors-{{neighbors}}-resolution-{{res}}.pdf', 
            cell_types = all_cell_types['CyTOF'], pn = ['positive']),
        ground_truth_umap_pdf = output + 'figures/CyTOF-Seurat-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}.pdf'
    script:
        'cell-type-assignment/Seurat.R'


rule CyTOF_LDA:
    input:
        training_rds = 'data/CyTOF/CyTOF-train.rds',
        annotation_rds = 'data/CyTOF/CyTOF-test.rds',
        labels = lambda wildcards: get_labels(wildcards.selection_procedure, 'CyTOF'),
    output:
        prediction = output + 'rare-subtype-benchmarking/CyTOF-{selection_procedure}-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-{annotator}-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-{set}-CyTOF-LDA-predictions-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/CyTOFLDA.R'
        
