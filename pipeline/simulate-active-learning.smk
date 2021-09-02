

active_learner = {
    'Create_labels': expand('data/{modality}/Active-Learning/Active-Learning-annotation-GroundTruth.tsv', modality = modalities),
    'entropy_list': expand('data/{modality}/Active-Learning/Active-Learning-entropy-GroundTruth.tsv', modality = modalities),
    'entropy_box': expand(output + '{modality}-Active-Learning-entropy-GroundTruth-boxplot.pdf', modality = modalities),
    'entropy_hist': expand(output + '{modality}-Active-Learning-entropy-GroundTruth-histogram.pdf', modality = modalities),
    'subset_al': expand('data/{modality}/Active-Learning/Active-Learning-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{subset_val}.tsv', 
        modality = modalities, subset_val = list(range(15)))
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
        box = output + '{modality}-Active-Learning-entropy-GroundTruth-boxplot.pdf',
        hist = output + '{modality}-Active-Learning-entropy-GroundTruth-histogram.pdf',
    script:
        'visualize/plot-entropies.R'


rule subset_training_data:
    input:
        assignment = 'data/{modality}/Active-Learning/Active-Learning-annotation-GroundTruth.tsv'
    output:
        split = 'data/{modality}/Active-Learning/Active-Learning-{modality}-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{subset_val}.tsv'
    script:
        'cell-type-assignment/subset-simulated-active-learner.R'