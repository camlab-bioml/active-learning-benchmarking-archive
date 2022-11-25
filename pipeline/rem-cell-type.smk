
rem_cell_type_reports = []
rctreport = [expand(output + 'results/rem_cell_type/Init-sel-{initial}-rem-{rem_cell_type}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}-seed-{s}.tsv',
    initial = initial_selections,
    rem_cell_type = original_cell_types[modality],
    modality = [modality],
    AL_alg = AL_methods,
    strat = ['highest_entropy', 'lowest_maxp'],
    num = [20, 50],
    s = train_test_seeds)
    for modality in original_cell_types.keys()]

for element in rctreport:
    rem_cell_type_reports.extend(element)

rem_cell_type = {
    #'reports': rem_cell_type_reports,
    #'summary': output + 'results/rem_cell_type/Rem-cell-types-combined-uncertanties.tsv',
    'plot': expand(output + 'figures/rem_cell_type/Selected_cells-Init-sel-{initial}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}.pdf',
        initial = initial_selections, modality = modalities, strat = ['highest_entropy', 'lowest_maxp'], AL_alg = AL_methods, num = [20, 50]),
    'probs': expand(output + 'figures/rem_cell_type_prob/Selected_cells-Init-sel-{initial}-{modality}-ALAlg-{AL_alg}-cell_num-{num}-rem_celltype-{rem_cell_type}-seed-{s}.pdf',
        initial = ['random'], modality = ['scRNASeq'], AL_alg = AL_methods, num = [20, 50], rem_cell_type = original_cell_types['scRNASeq'], s = [0])
}

rule rem_cell_type:
    input:
        markers = 'markers/{modality}.yml',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    resources:
        mem_mb=10000
    output:
        tsv = output + 'results/rem_cell_type/Init-sel-{initial}-rem-{rem_cell_type}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}-seed-{s}.tsv'
    script:
        "rem-cell-type-from-training/rem-cell-type.R"

rule combine_uncertainties:
    input:
        tsv = rem_cell_type_reports
    output:
        tsv = output + 'results/rem_cell_type/Rem-cell-types-combined-uncertanties.tsv'
    script:
        'rem-cell-type-from-training/combine-uncertainties.R'

def get_uncertainties_by_modality(mod, mod_cts = original_cell_types, seeds = train_test_seeds):
    cts = mod_cts[mod]

    f = expand(output + 'results/rem_cell_type/Init-sel-{{initial}}-rem-{rem_cell_type}-{{modality}}-{{strat}}-ALAlg-{{AL_alg}}-cell_num-{{num}}-seed-{s}.tsv',
        rem_cell_type = cts, s = seeds)
    
    return f


rule plot_proportion_selected:
    input:
        tsvs = lambda wildcards: get_uncertainties_by_modality(wildcards.modality)
    output:
        plot = output + 'figures/rem_cell_type/Selected_cells-Init-sel-{initial}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}.pdf'
    script:
        'rem-cell-type-from-training/plot-proportion-selected.R'

rule understand_probabilities:
    input:
        markers = 'markers/{modality}.yml',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    output:
        pdf = output + 'figures/rem_cell_type_prob/Selected_cells-Init-sel-{initial}-{modality}-ALAlg-{AL_alg}-cell_num-{num}-rem_celltype-{rem_cell_type}-seed-{s}.pdf'
    script:
        'rem-cell-type-from-training/understand-probs.R'