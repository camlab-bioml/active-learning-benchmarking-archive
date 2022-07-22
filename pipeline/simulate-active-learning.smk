

active_learner = {
    'Create_labels': expand('data/{modality}/Active-Learning/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-full.tsv', 
        modality = modalities, strat = selection_expansion_dict['Active-Learning']['strategy'], rand = selection_expansion_dict['Active-Learning']['random_selection'], corrupt = selection_expansion_dict['Active-Learning']['corruption']),
    #'entropy_list': expand('data/{modality}/Active-Learning/Active-Learning-entropy-GroundTruth.tsv', modality = modalities),
    #'entropy_box': expand(output + 'figures/AL-benchmarking/{modality}-Active-Learning-entropy-GroundTruth-boxplot.pdf', modality = modalities),
    #'entropy_hist': expand(output + 'figures/AL-benchmarking/{modality}-Active-Learning-entropy-GroundTruth-histogram.pdf', modality = modalities),
    'subset_al': expand('data/{modality}/Active-Learning/AL-batches-subset/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{subset_val}_cells.tsv', 
        modality = modalities, strat = selection_expansion_dict['Active-Learning']['strategy'], rand = selection_expansion_dict['Active-Learning']['random_selection'], corrupt = selection_expansion_dict['Active-Learning']['corruption'], subset_val = cell_numbers),
    #'removed_marker_AL': expand(output + 'AL_with_removed_input/scRNASeq-marker-{removed}-seed-{seed}.tsv', removed = ['removed', 'kept'], seed = list(range(10)))
}


rule Create_active_learning_ground_truth:
    input:
        markers = 'markers/{modality}.yml',
        expression = 'data/{modality}/{modality}-train.rds'
    output:
        assignments = 'data/{modality}/Active-Learning/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-full.tsv',
        entropy = 'data/{modality}/Active-Learning/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-entropy-GroundTruth.tsv'
    script:
        'cell-type-assignment/simulate-active-learner.R'


rule visualize_entropies:
    input:
        entropy = 'data/{modality}/Active-Learning/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-entropy-GroundTruth.tsv'
    output:
        box = output + 'figures/AL-benchmarking/{modality}-Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-entropy-GroundTruth-boxplot.pdf',
        hist = output + 'figures/AL-benchmarking/{modality}-Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-entropy-GroundTruth-histogram.pdf',
    script:
        'visualize/plot-entropies.R'


rule subset_training_data:
    # This rule subsets the training data by selecting batches of ranked cells
    input:
        assignment = 'data/{modality}/Active-Learning/Active-Learning-{strat}-{modality}-rand_sel-{rand}-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-full.tsv'
    output:
        split = 'data/{modality}/Active-Learning/AL-batches-subset/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{subset_val}_batches.tsv'
    script:
        'cell-type-assignment/subset-simulated-active-learner.R'

rule create_AL_training_batches:
    # This rule subsets the training data by selecting the first n cells
    input:
        assignment = 'data/{modality}/Active-Learning/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-full.tsv'
    output:
        split = 'data/{modality}/Active-Learning/AL-batches-subset/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{subset_val}_cells.tsv'
    script:
        'cell-type-assignment/subset-simulated-active-learner-n-cells.R'


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