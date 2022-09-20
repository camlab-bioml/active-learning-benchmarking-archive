

active_learner = {
    # 'Create_labels_entropy': expand(output + 'data/{modality}/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv', 
    #     modality = modalities, 
    #     AL_type = ['Active-Learning_entropy'], 
    #     initial = selection_expansion_dict['Active-Learning_entropy']['initial_sel'],
    #     strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], 
    #     AL_alg = selection_expansion_dict['Active-Learning_entropy']['AL_alg'],
    #     rand = selection_expansion_dict['Active-Learning_entropy']['random_selection'], 
    #     corrupt = selection_expansion_dict['Active-Learning_entropy']['corruption'],
    #     s = train_test_seeds),
    # 'Create_labels_maxp': expand(output + 'data/{modality}/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv', 
    #     modality = modalities, 
    #     AL_type = ['Active-Learning_maxp'], 
    #     initial = selection_expansion_dict['Active-Learning_maxp']['initial_sel'],
    #     strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], 
    #     AL_alg = selection_expansion_dict['Active-Learning_maxp']['AL_alg'],
    #     rand = selection_expansion_dict['Active-Learning_maxp']['random_selection'], 
    #     corrupt = selection_expansion_dict['Active-Learning_maxp']['corruption'],
    #     s = train_test_seeds),

    'maxp_box_rand_0': expand(output + 'figures/uncertainty-{modality}-Init-{initial}-strat-{AL_type}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.pdf',
        modality = modalities, AL_type = ['Active-Learning_maxp', 'Active-Learning_entropy'], initial = initial_selections, 
        AL_alg = AL_methods,
        rand = 0, 
        corrupt = corruption_percentages,
        s = train_test_seeds),
    'maxp_box_cor_0': expand(output + 'figures/uncertainty-{modality}-Init-{initial}-strat-{AL_type}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.pdf',
        modality = modalities, AL_type = ['Active-Learning_maxp', 'Active-Learning_entropy'], initial = initial_selections, 
        AL_alg = AL_methods,
        rand = random_percentages, 
        corrupt = 0,
        s = train_test_seeds),

    'subset_al_entropy_rand_0': expand(output + 'data/{modality}/{AL_type}/AL-batches-subset/Init-{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-{subset_val}_cells.tsv', 
        modality = modalities, initial = selection_expansion_dict['Active-Learning_entropy']['initial_sel'], AL_type = ['Active-Learning_entropy'], strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], 
        AL_alg = selection_expansion_dict['Active-Learning_entropy']['AL_alg'],
        rand = 0, corrupt = selection_expansion_dict['Active-Learning_entropy']['corruption'], 
        subset_val = cell_numbers, s = train_test_seeds),

    'subset_al_entropy_corr_0': expand(output + 'data/{modality}/{AL_type}/AL-batches-subset/Init-{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-{subset_val}_cells.tsv', 
        modality = modalities, initial = selection_expansion_dict['Active-Learning_entropy']['initial_sel'], AL_type = ['Active-Learning_entropy'], strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], 
        AL_alg = selection_expansion_dict['Active-Learning_entropy']['AL_alg'],
        rand = selection_expansion_dict['Active-Learning_entropy']['random_selection'], corrupt = 0, 
        subset_val = cell_numbers, s = train_test_seeds),

    'subset_al_maxp_corr_0': expand(output + 'data/{modality}/{AL_type}/AL-batches-subset/Init-{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-{subset_val}_cells.tsv', 
        modality = modalities, initial = selection_expansion_dict['Active-Learning_maxp']['initial_sel'], AL_type = ['Active-Learning_maxp'], strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], 
        AL_alg = selection_expansion_dict['Active-Learning_maxp']['AL_alg'],
        rand = selection_expansion_dict['Active-Learning_maxp']['random_selection'], corrupt = 0, 
        subset_val = cell_numbers, s = train_test_seeds),
    'subset_al_maxp_rand_0': expand(output + 'data/{modality}/{AL_type}/AL-batches-subset/Init-{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-{subset_val}_cells.tsv', 
        modality = modalities, initial = selection_expansion_dict['Active-Learning_maxp']['initial_sel'], AL_type = ['Active-Learning_maxp'], strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], 
        AL_alg = selection_expansion_dict['Active-Learning_maxp']['AL_alg'],
        rand = 0, corrupt = selection_expansion_dict['Active-Learning_maxp']['corruption'], 
        subset_val = cell_numbers, s = train_test_seeds),
    #'removed_marker_AL': expand(output + 'AL_with_removed_input/scRNASeq-marker-{removed}-seed-{seed}.tsv', removed = ['removed', 'kept'], seed = list(range(10)))
}

## CHECKED
rule Create_active_learning_ground_truth:
    input:
        markers = 'markers/{modality}.yml',
        expression = 'data/{modality}/{modality}-train-seed-{s}.rds'
    params:
        max_cell_num = max(cell_numbers)
    output:
        assignments = output + 'data/{modality}/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv',
        entropy = output + 'data/{modality}/uncertainties/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv'
    script:
        'cell-type-assignment/simulate-active-learner.R'


def get_entropy_files(AL_type):
    if AL_type == "Active-Learning_entropy":
        f = expand(output + 'data/{{modality}}/uncertainties/Active-Learning_entropy/Init-{{initial}}-strat-{strat}-al-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-seed-{{s}}.tsv',
            strat = selection_expansion_dict['Active-Learning_entropy']['strategy'])
    elif AL_type == 'Active-Learning_maxp':
        f = expand(output + 'data/{{modality}}/uncertainties/Active-Learning_maxp/Init-{{initial}}-strat-{strat}-al-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-seed-{{s}}.tsv',
            strat = selection_expansion_dict['Active-Learning_maxp']['strategy'])
    return f

## CHECKED
rule visualize_entropies:
    input:
        f = lambda wildcards: get_entropy_files(wildcards.AL_type)
    output:
        box = output + 'figures/uncertainty-{modality}-Init-{initial}-strat-{AL_type}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.pdf'
    script:
        'visualize/plot-entropies.R'


# rule subset_training_data:
#     # This rule subsets the training data by selecting batches of ranked cells
#     input:
#         assignment = 'data/{modality}/{initial}_{AL_type}/{AL_type}-{strat}-{modality}-rand_sel-{rand}-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-full.tsv'
#     output:
#         split = 'data/{modality}/{initial}_{AL_type}/AL-batches-subset/{AL_type}-{strat}-rand_sel-{rand}-corr-{corrupt}-{modality}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{subset_val}_batches.tsv'
#     script:
#         'cell-type-assignment/subset-simulated-active-learner.R'

## CHECKED ###
rule create_AL_training_batches:
    # This rule subsets the training data by selecting the first n cells
    input:
        assignment = output + 'data/{modality}/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv'
    output:
        split = output + 'data/{modality}/{AL_type}/AL-batches-subset/Init-{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-{subset_val}_cells.tsv'
    script:
        'cell-type-assignment/subset-simulated-active-learner-n-cells.R'


# rule AL_with_removed_input:
#     input:
#         markers = 'markers/scRNASeq.yml',
#         expression = 'data/scRNASeq/scRNASeq-train.rds'
#     output:
#         tsv = output + 'AL_with_removed_input/scRNASeq-marker-removed-seed-{seed}.tsv'
#     script:
#         'benchmarking/remove-cellType-from-ranked-cells.R'

# rule AL_with_kept_markers:
#     input:
#         markers = 'markers/scRNASeq.yml',
#         expression = 'data/scRNASeq/scRNASeq-train.rds'
#     output:
#         tsv = output + 'AL_with_removed_input/scRNASeq-marker-kept-seed-{seed}.tsv'
#     script:
#         'benchmarking/keep-cellType-from-ranked-cells.R'