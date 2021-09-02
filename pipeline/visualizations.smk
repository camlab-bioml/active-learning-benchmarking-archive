

viz = {
    'AL_iterations': expand(output + 'figures/AL-benchmarking/Active-Learning-accuracy-by-training-iteration-{modality}.pdf', modality = ['scRNASeq'])
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