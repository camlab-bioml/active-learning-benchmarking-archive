

viz = {
    #'AL_iterations': expand(output + 'figures/AL-benchmarking/Active-Learning-accuracy-by-training-iteration-{modality}.pdf', modality = ['scRNASeq'])
    'benchmark': expand(output + 'reports/overall-{modality}-benchmarking-cells.html', modality = modalities),
    'selected_cells': expand(output + 'reports/vizualize-selected-{modality}-cells.html', modality = modalities),
    'viz_corruption': expand(output + 'reports/vizualize-corruption-{modality}-cells.html', modality = modalities),
    'gt_predicted_alluv': expand(output + 'reports/predicted-ground-truth-alluvials-{modality}.html', modality = modalities)
}

def create_AL_iterations_assignments(mod):
    if(mod == 'scRNASeq'):
        methods = evaluation_methods_dict['scRNASeq']
    elif(mod == 'CyTOF'):
        methods = evaluation_methods_dict['CyTOF']

    return(expand(output + 'rare-subtype-benchmarking/{{modality}}-Active-Learning-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{set}-{method}-predictions.tsv',
           set = list(range(training_set_AL_simulation)), method = methods))

rule AL_iterations_benchmark:
    input:
        sce = "data/{modality}/{modality}-test.rds",
        assignments = lambda wildcards: create_AL_iterations_assignments(wildcards.modality)
    output:
        pdf = output + 'figures/AL-benchmarking/Active-Learning-accuracy-by-training-iteration-{modality}.pdf'
    script:
       'visualize/AL-iterations-benchmark.R' 


rule visualize_selected_cells:
    input:
        markers = 'markers/{modality}.yml',
        random_files = expand('data/{{modality}}/random/random-NA-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{set}-{cell_num}_cells.tsv',
            rand = selection_expansion_dict['random']['random_selection'], corrupt = selection_expansion_dict['random']['corruption'], set = random_sets, cell_num = cell_numbers),
        seurat_files = expand('data/{{modality}}/Seurat-clustering/Seurat-clustering-NA-rand_sel-NA-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-{neighbors}-resolution-{res}-iterations_set-NA-{cell_num}_cells.tsv',
            corrupt = selection_expansion_dict['Seurat-clustering']['corruption'], neighbors = selection_expansion_dict['Seurat-clustering']['neighbors'], res = selection_expansion_dict['Seurat-clustering']['res'], cell_num = cell_numbers),
        AL = expand('data/{{modality}}/Active-Learning/AL-batches-subset/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{cell_num}_cells.tsv',
            strat = selection_expansion_dict['Active-Learning']['strategy'], rand = selection_expansion_dict['Active-Learning']['random_selection'], corrupt = selection_expansion_dict['Active-Learning']['corruption'], cell_num = cell_numbers),
        sce = 'data/{modality}/{modality}-train.rds'
    params:
        output_dir = output + 'reports/',
        random_files = 'data/{modality}/random/',
        seurat_files = 'data/{modality}/Seurat-clustering/',
        AL_files = 'data/{modality}/Active-Learning/AL-batches-subset/',
    output:
        html = output + 'reports/vizualize-selected-{modality}-cells.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/selected-cells-by-method.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', markers = '{input.markers}', random_files_path = '{params.random_files}', seurat_files_path = '{params.seurat_files}', AL_files_path = '{params.AL_files}', modality = '{wildcards.modality}'))\" "

rule visualize_corruption:
    input:
        AL = expand('data/{{modality}}/Active-Learning/AL-batches-subset/Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-{{modality}}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{cell_num}_cells.tsv',
            strat = selection_expansion_dict['Active-Learning']['strategy'], rand = selection_expansion_dict['Active-Learning']['random_selection'], corrupt = selection_expansion_dict['Active-Learning']['corruption'], cell_num = cell_numbers),
        sce = 'data/{modality}/{modality}-train.rds'
    params:
        output_dir = output + 'reports/',
        AL_files = 'data/{modality}/Active-Learning/AL-batches-subset/',
    output:
        html = output + 'reports/vizualize-corruption-{modality}-cells.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/visualize/visualize-corruption.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', AL_files_path = '{params.AL_files}', modality = '{wildcards.modality}'))\" "

def get_predictions(mod, selection):
    if mod == 'scRNASeq':
        if selection == "AL":
            f = expand(output + 'rare-subtype-benchmarking/scRNASeq-Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{method}-predictions-{cell_num}-cells.tsv', 
                strat = selection_expansion_dict['Active-Learning']['strategy'], rand = selection_expansion_dict['Active-Learning']['random_selection'], corrupt = selection_expansion_dict['Active-Learning']['corruption'], method = evaluation_methods_dict['scRNASeq'], cell_num = cell_numbers)
        elif selection == 'Seurat':
            f = expand(output + 'rare-subtype-benchmarking/scRNASeq-Seurat-clustering-NA-rand_sel-NA-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-{knn}-resolution-{res}-iterations_set-NA-{method}-predictions-{cell_num}-cells.tsv',
                corrupt = selection_expansion_dict['Seurat-clustering']['corruption'], knn = Seurat_neighbors, res = Seurat_resolution, method = evaluation_methods_dict['scRNASeq'], cell_num = cell_numbers)
        elif selection == 'random':
            f =  expand(output + 'rare-subtype-benchmarking/scRNASeq-random-NA-rand_sel-NA-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{set}-{method}-predictions-{cell_num}-cells.tsv', 
                corrupt = selection_expansion_dict['random']['corruption'], set = random_sets, method = evaluation_methods_dict['scRNASeq'], cell_num = cell_numbers)
    elif mod == 'CyTOF':
        if selection == "AL":
            f = expand(output + 'rare-subtype-benchmarking/CyTOF-Active-Learning-{strat}-rand_sel-{rand}-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA-{method}-predictions-{cell_num}-cells.tsv',
                strat = selection_expansion_dict['Active-Learning']['strategy'], rand = selection_expansion_dict['Active-Learning']['random_selection'], corrupt = selection_expansion_dict['Active-Learning']['corruption'], method = evaluation_methods_dict['CyTOF'], cell_num = cell_numbers)
        elif selection == 'Seurat':
            f = expand(output + 'rare-subtype-benchmarking/CyTOF-Seurat-clustering-NA-rand_sel-NA-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-{knn}-resolution-{res}-iterations_set-NA-{method}-predictions-{cell_num}-cells.tsv',
                corrupt = selection_expansion_dict['Seurat-clustering']['corruption'], knn = Seurat_neighbors, res = Seurat_resolution, method = evaluation_methods_dict['CyTOF'], cell_num = cell_numbers)
        elif selection == 'random':
            f =  expand(output + 'rare-subtype-benchmarking/CyTOF-random-NA-rand_sel-NA-corr-{corrupt}-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{set}-{method}-predictions-{cell_num}-cells.tsv', 
                corrupt = selection_expansion_dict['random']['corruption'], set = random_sets, method = evaluation_methods_dict['CyTOF'], cell_num = cell_numbers)
    
    return (f)

rule overall_benchmark:
    input:
        sce = 'data/{modality}/{modality}-test.rds',
        active_learning = lambda wildcards: get_predictions(wildcards.modality, "AL"),
        seurat = lambda wildcards: get_predictions(wildcards.modality, "Seurat"),
        random = lambda wildcards: get_predictions(wildcards.modality, "random")
    params:
        input_dir = output + 'rare-subtype-benchmarking',
        AL_pattern = "{modality}-Active-Learning-*annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA*-cells*",
        seurat_pattern = "{modality}-Seurat*-cells*",
        random_pattern = "{modality}-random*-cells*",
        output_dir = output + 'reports/'
    output:
        html = output + 'reports/overall-{modality}-benchmarking-cells.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/benchmarking.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', input_dir = '{params.input_dir}', AL_pattern = '{params.AL_pattern}', "
        "seurat_pattern = '{params.seurat_pattern}', random_pattern = '{params.random_pattern}'))\" "


rule compare_predicted_to_ground_truth_alluvial:
    input:
        sce = 'data/{modality}/{modality}-test.rds',
        active_learning = lambda wildcards: get_predictions(wildcards.modality, "AL"),
        seurat = lambda wildcards: get_predictions(wildcards.modality, "Seurat"),
        random = lambda wildcards: get_predictions(wildcards.modality, "random")
    params:
        input_dir = output + 'rare-subtype-benchmarking/',
        AL_pattern = "{modality}-Active-Learning-*annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-NA*-cells*",
        seurat_pattern = "{modality}-Seurat*-cells*",
        random_pattern = "{modality}-random*-cells*",
        output_dir = output + 'reports/'
    output:
        html = output + 'reports/predicted-ground-truth-alluvials-{modality}.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/visualize/predicted-ground-truth-alluvials.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', input_dir = '{params.input_dir}', AL_pattern = '{params.AL_pattern}', "
        "seurat_pattern = '{params.seurat_pattern}', random_pattern = '{params.random_pattern}'))\" "
