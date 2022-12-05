

similarity = {
    'heatmaps': expand(output + 'figures/cell-type-similarity/similarity-heatmap-{modality}-seed-{s}.pdf', modality = modalities, s = train_test_seeds),
    'tsvs': expand(output + 'results/cell-type-similarity/similarity-{modality}-seed-{s}.tsv', modality = modalities, s = train_test_seeds),
    'mean': expand(output + 'results/cell-type-similarity/average-similarity-{modality}.tsv', modality = modalities)
}


rule calculate_cosine_dist:
    input:
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    output:
        heatmap = output + 'figures/cell-type-similarity/similarity-heatmap-{modality}-seed-{s}.pdf',
        tsv = output + 'results/cell-type-similarity/similarity-{modality}-seed-{s}.tsv'
    script:
        'cell-type-similarity/calculate-cosine-distance.R'

rule summarize_cosine_dist_across_seeds:
    input:
        tsvs = expand(output + 'results/cell-type-similarity/similarity-{modality}-seed-{s}.tsv',
            modality = modalities, s = train_test_seeds)
    output:
        avg_dist_sc = output + 'results/cell-type-similarity/average-similarity-scRNASeq.tsv',
        avg_dist_sn = output + 'results/cell-type-similarity/average-similarity-snRNASeq.tsv',
        avg_dist_cy = output + 'results/cell-type-similarity/average-similarity-CyTOF.tsv'
    script:
        'cell-type-similarity/summarize-cosine-distances.R'