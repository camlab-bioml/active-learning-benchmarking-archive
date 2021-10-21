

viz = {
    #'AL_iterations': expand(output + 'figures/AL-benchmarking/Active-Learning-accuracy-by-training-iteration-{modality}.pdf', modality = ['scRNASeq'])
    'benchmark': expand(output + 'reports/overall-{modality}-benchmarking.html', modality = ['scRNASeq']),
    'selected_cells': expand(output + 'reports/vizualize-selected-{modality}-cells.html', modality = ['scRNASeq']),
    'gt_predicted_alluv': output + 'reports/predicted-ground-truth-alluvials.html'
}


rule AL_iterations_benchmark:
    input:
        sce = "data/{modality}/{modality}-test.rds",
        assignments = expand(output + 'rare-subtype-benchmarking/{{modality}}-Active-Learning-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{set}-{method}-predictions.tsv',
                                                                     #scRNASeq-Active-Learning-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-0-scmap-cluster-predictions
                               set = list(range(training_set_AL_simulation)), method = scRNASeq_methods)
    output:
        pdf = output + 'figures/AL-benchmarking/Active-Learning-accuracy-by-training-iteration-{modality}.pdf'
    script:
       'visualize/AL-iterations-benchmark.R' 


rule visualize_selected_cells:
    input:
        markers = 'markers/scRNASeq.yml',
        random_files = expand('data/scRNASeq/random/random-scRNASeq-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{set}.tsv', set = random_sets),
        seurat_files = expand('data/scRNASeq/Seurat-clustering/Seurat-clustering-scRNASeq-annotator-GroundTruth-knn_neighbors-{knn}-resolution-{res}-iterations_set-NA.tsv', knn = Seurat_neighbors, res = Seurat_resolution),
        AL_file = 'data/scRNASeq/Active-Learning/Active-Learning-annotation-GroundTruth.tsv',
        sce = 'data/scRNASeq/scRNASeq-train.rds'
    params:
        output_dir = output + 'reports/'
    output:
        html = output + 'reports/vizualize-selected-scRNASeq-cells.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/selected-cells-by-method.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', markers = '{input.markers}', random_files = '{input.random_files}', seurat_files = '{input.seurat_files}', AL_file = '{input.AL_file}'))\" "


rule overall_benchmark:
    input:
        sce = 'data/scRNASeq/scRNASeq-test.rds',
        active_learning = expand(output + 'rare-subtype-benchmarking/scRNASeq-Active-Learning-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-14-{method}-predictions.tsv', method = scRNASeq_methods),
        seurat = expand(output + 'rare-subtype-benchmarking/scRNASeq-Seurat-clustering-annotator-GroundTruth-knn_neighbors-{knn}-resolution-{res}-iterations_set-NA-{method}-predictions.tsv', knn = Seurat_neighbors, res = Seurat_resolution, method = scRNASeq_methods),
        random = expand(output + 'rare-subtype-benchmarking/scRNASeq-random-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{set}-{method}-predictions.tsv', set = random_sets, method = scRNASeq_methods)
    params:
        input_dir = output + 'rare-subtype-benchmarking/',
        AL_pattern = "scRNASeq-Active-Learning-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-14",
        seurat_pattern = "scRNASeq-Seurat",
        random_pattern = "scRNASeq-random",
        output_dir = output + 'reports/'
    output:
        html = output + 'reports/overall-{modality}-benchmarking.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/benchmarking.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', input_dir = '{params.input_dir}', AL_pattern = '{params.AL_pattern}', "
        "seurat_pattern = '{params.seurat_pattern}', random_pattern = '{params.random_pattern}'))\" "


rule compare_predicted_to_ground_truth_alluvial:
    input:
        sce = 'data/scRNASeq/scRNASeq-test.rds',
        active_learning = expand(output + 'rare-subtype-benchmarking/scRNASeq-Active-Learning-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-14-{method}-predictions.tsv', method = scRNASeq_methods),
        seurat = expand(output + 'rare-subtype-benchmarking/scRNASeq-Seurat-clustering-annotator-GroundTruth-knn_neighbors-{knn}-resolution-{res}-iterations_set-NA-{method}-predictions.tsv', knn = Seurat_neighbors, res = Seurat_resolution, method = scRNASeq_methods),
        random = expand(output + 'rare-subtype-benchmarking/scRNASeq-random-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-{set}-{method}-predictions.tsv', set = random_sets, method = scRNASeq_methods)
    params:
        input_dir = output + 'rare-subtype-benchmarking/',
        AL_pattern = "scRNASeq-Active-Learning-annotator-GroundTruth-knn_neighbors-NA-resolution-NA-iterations_set-14",
        seurat_pattern = "scRNASeq-Seurat",
        random_pattern = "scRNASeq-random",
        output_dir = output + 'reports/'
    output:
        html = output + 'reports/predicted-ground-truth-alluvials.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/visualize/predicted-ground-truth-alluvials.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', input_dir = '{params.input_dir}', AL_pattern = '{params.AL_pattern}', "
        "seurat_pattern = '{params.seurat_pattern}', random_pattern = '{params.random_pattern}'))\" "
