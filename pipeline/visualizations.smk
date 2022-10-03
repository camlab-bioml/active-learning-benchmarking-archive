

viz = {
    # 'benchmark': expand(output + 'reports/overall-{modality}-{AL_alg}-benchmarking-cells.html', modality = modalities, AL_alg = AL_methods),
    'cross_cohort': output + 'reports/cross-cohort-eval.html',
    #'selected_cells': expand(output + 'reports/vizualize-selected-{modality}-cells.html', modality = modalities),
    #'viz_corruption': expand(output + 'reports/vizualize-corruption-{modality}-cells.html', modality = modalities),
    #'gt_predicted_alluv': expand(output + 'reports/predicted-ground-truth-alluvials-{modality}.html', modality = modalities)
}

def create_AL_iterations_assignments(mod):
    if(mod == 'scRNASeq'):
        methods = evaluation_methods_dict['scRNASeq']
    elif(mod == 'CyTOF'):
        methods = evaluation_methods_dict['CyTOF']

    return(expand(output + 'rare-subtype-benchmarking/{{modality}}-Active-Learning-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{set}-{method}-predictions.tsv',
           set = list(range(training_set_AL_simulation)), method = methods))


rule visualize_selected_cells:
    input:
        markers = 'markers/{modality}.yml',
        random_files = expand(output + 'data/{{modality}}/random/random-NA-ALAlg-NA-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{set}-{cell_num}_cells.tsv',
            rand = selection_expansion_dict['random']['random_selection'], corrupt = selection_expansion_dict['random']['corruption'], set = random_sets, cell_num = cell_numbers),
        seurat_files = expand(output + 'data/{{modality}}/Seurat-clustering/Seurat-clustering-NA-ALAlg-NA-rand_sel-NA-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-NA-{cell_num}_cells.tsv',
            corrupt = selection_expansion_dict['Seurat-clustering']['corruption'], neighbors = selection_expansion_dict['Seurat-clustering']['neighbors'], res = selection_expansion_dict['Seurat-clustering']['res'], cell_num = cell_numbers),
        AL_entropy = expand(output + 'data/{{modality}}/Active-Learning_entropy/AL-batches-subset/Active-Learning_entropy-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{cell_num}_cells.tsv',
            strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], AL_alg = AL_methods, rand = selection_expansion_dict['Active-Learning_entropy']['random_selection'], corrupt = selection_expansion_dict['Active-Learning_entropy']['corruption'], cell_num = cell_numbers),
        AL_maxp = expand(output + 'data/{{modality}}/Active-Learning_maxp/AL-batches-subset/Active-Learning_maxp-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{cell_num}_cells.tsv',
            strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], AL_alg = AL_methods, rand = selection_expansion_dict['Active-Learning_maxp']['random_selection'], corrupt = selection_expansion_dict['Active-Learning_maxp']['corruption'], cell_num = cell_numbers),
        sce = 'data/{modality}/{modality}-train.rds'
    params:
        output_dir = output + 'reports/',
        random_files = output + 'data/{modality}/random/',
        seurat_files = output + 'data/{modality}/Seurat-clustering/',
        AL_files_entropy = output + 'data/{modality}/Active-Learning_entropy/AL-batches-subset/',
        AL_files_maxp = output + 'data/{modality}/Active-Learning_maxp/AL-batches-subset/'
    output:
        html = output + 'reports/vizualize-selected-{modality}-cells.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/selected-cells-by-method.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', markers = '{input.markers}', random_files_path = '{params.random_files}', seurat_files_path = '{params.seurat_files}', "
        "AL_entropy_path = '{params.AL_files_entropy}', AL_maxp_path = '{params.AL_files_maxp}', modality = '{wildcards.modality}'))\" "

rule visualize_corruption:
    input:
        AL_entropy = expand(output + 'data/{{modality}}/Active-Learning_entropy/AL-batches-subset/Active-Learning_entropy-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{cell_num}_cells.tsv',
            AL_alg = AL_methods, strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], rand = selection_expansion_dict['Active-Learning_entropy']['random_selection'], corrupt = selection_expansion_dict['Active-Learning_entropy']['corruption'], cell_num = cell_numbers),
        AL_maxp = expand(output + 'data/{{modality}}/Active-Learning_maxp/AL-batches-subset/Active-Learning_maxp-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{cell_num}_cells.tsv',
            AL_alg = AL_methods, strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], rand = selection_expansion_dict['Active-Learning_maxp']['random_selection'], corrupt = selection_expansion_dict['Active-Learning_maxp']['corruption'], cell_num = cell_numbers),
        sce = 'data/{modality}/{modality}-train.rds'
    params:
        output_dir = output + 'reports/',
        AL_files_entropy = output + 'data/{modality}/Active-Learning_entropy/AL-batches-subset/',
        AL_files_maxp = output + 'data/{modality}/Active-Learning_maxp/AL-batches-subset/'
    output:
        html = output + 'reports/vizualize-corruption-{modality}-cells.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/visualize/visualize-corruption.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', AL_entropy_path = '{params.AL_files_entropy}', AL_maxp_path = '{params.AL_files_maxp}', modality = '{wildcards.modality}'))\" "

def get_predictions(mod):
    if mod == 'scRNASeq':
        f = cell_type_predictions['scRNASeq']
    elif mod == 'snRNASeq':
        f = cell_type_predictions['snRNASeq']
    elif mod == 'CyTOF':
        f = cell_type_predictions['CyTOF']
    
    return (f)

rule overall_benchmark:
    input:
        sce = 'data/{modality}/{modality}-full.rds',
        predictions = lambda wildcards: get_predictions(wildcards.modality)
    resources:
        mem_mb=50000
    log:
        output + 'logs/benchmark-predictive-labeling-{modality}.log'
    output:
        acc = output + 'results/overall-{modality}-benchmarking-accuracies.tsv'
    script:
        'benchmarking/save-acc-overall-benchmarking.R'

rule visualize_overall_benchmark:
    input:
        acc = output + 'results/overall-{modality}-benchmarking-accuracies.tsv'
    output:
        html = output + 'reports/overall-{modality}-{AL_alg}-benchmarking-cells.html'
    params:
        output_dir = output + 'reports/'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/benchmarking.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(acc = '{input.acc}', al_alg = '{wildcards.AL_alg}'))\" "

rule cross_cohort_overview:
    input:
        accs = expand(output + 'results/overall-{modality}-benchmarking-accuracies.tsv', modality = modalities)
    output:
        html = output + 'reports/cross-cohort-eval.html'
    params:
        output_dir = output + 'reports/'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/Cross-cohort-eval.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(accs = '{input.accs}'))\" "

rule compare_predicted_to_ground_truth_alluvial:
    input:
        sce = 'data/{modality}/{modality}-test.rds',
        active_learning_entropy = lambda wildcards: get_predictions(wildcards.modality, "AL_entropy"),
        active_learning_maxp = lambda wildcards: get_predictions(wildcards.modality, "AL_maxp"),
        seurat = lambda wildcards: get_predictions(wildcards.modality, "Seurat"),
        random = lambda wildcards: get_predictions(wildcards.modality, "random")
    params:
        input_dir = output + 'rare-subtype-benchmarking/',
        AL_pattern_entropy = "{modality}-sel-Active-Learning_entropy-*annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA*-cells*",
        AL_pattern_maxp = "{modality}-sel-Active-Learning_maxp-*annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA*-cells*",
        seurat_pattern = "{modality}-sel-Seurat*-cells*",
        random_pattern = "{modality}-sel-random*-cells*",
        output_dir = output + 'reports/'
    output:
        html = output + 'reports/predicted-ground-truth-alluvials-{modality}.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/visualize/predicted-ground-truth-alluvials.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', input_dir = '{params.input_dir}', AL_pattern_entropy = '{params.AL_pattern_entropy}', "
        "AL_pattern_maxp = '{params.AL_pattern_maxp}', seurat_pattern = '{params.seurat_pattern}', random_pattern = '{params.random_pattern}', modality = '{wildcards.modality}'))\" "
