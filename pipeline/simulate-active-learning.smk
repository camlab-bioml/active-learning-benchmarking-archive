

active_learner = {
    #'Create_labels': expand('data/{modality}/Active-Learning/Active-Learning-annotation-GroundTruth.tsv', modality = modalities),
    #'entropy_list': expand('data/{modality}/Active-Learning/Active-Learning-entropy-GroundTruth.tsv', modality = modalities),
    #'entropy_box': expand(output + 'figures/AL-benchmarking/{modality}-Active-Learning-entropy-GroundTruth-boxplot.pdf', modality = modalities),
    #'entropy_hist': expand(output + 'figures/AL-benchmarking/{modality}-Active-Learning-entropy-GroundTruth-histogram.pdf', modality = modalities),
    #'subset_al': expand('data/{modality}/Active-Learning/Active-Learning-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{subset_val}.tsv', 
    #    modality = modalities, subset_val = list(range(15))),
    'removed_marker_AL': expand(output + 'AL_with_removed_input/scRNASeq-marker-{removed}-seed-{seed}.tsv', removed = ['removed', 'kept'], seed = list(range(10)))
}


rule Create_active_learning_ground_truth:
    input:
        markers = 'markers/{modality}.yml',
        expression = 'data/{modality}/{modality}-train.rds'
    output:
        assignments = 'data/{modality}/Active-Learning/Active-Learning-annotation-GroundTruth.tsv',
        entropy = 'data/{modality}/Active-Learning/Active-Learning-entropy-GroundTruth.tsv'
    script:
        'cell-type-assignment/simulate-active-learner.R'


rule visualize_entropies:
    input:
        entropy = 'data/{modality}/Active-Learning/Active-Learning-entropy-GroundTruth.tsv'
    output:
        box = output + 'figures/AL-benchmarking/{modality}-Active-Learning-entropy-GroundTruth-boxplot.pdf',
        hist = output + 'figures/AL-benchmarking/{modality}-Active-Learning-entropy-GroundTruth-histogram.pdf',
    script:
        'visualize/plot-entropies.R'


rule subset_training_data:
    input:
        assignment = 'data/{modality}/Active-Learning/Active-Learning-annotation-GroundTruth.tsv'
    output:
        split = 'data/{modality}/Active-Learning/Active-Learning-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{subset_val}.tsv'
    script:
        'cell-type-assignment/subset-simulated-active-learner.R'


rule AL_with_removed_input:
    input:
        markers = 'markers/scRNASeq.yml',
        expression = 'data/scRNASeq/scRNASeq-train.rds'
    output:
        tsv = output + 'AL_with_removed_input/scRNASeq-marker-removed-seed-{seed}.tsv'
    script:
        'benchmarking/remove-cellType-from-ranked-cells.R'

rule AL_with_kept_markers:
    input:
        markers = 'markers/scRNASeq.yml',
        expression = 'data/scRNASeq/scRNASeq-train.rds'
    output:
        tsv = output + 'AL_with_removed_input/scRNASeq-marker-kept-seed-{seed}.tsv'
    script:
        'benchmarking/keep-cellType-from-ranked-cells.R'