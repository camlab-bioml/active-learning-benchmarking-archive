def create_AL_iterations_assignments(mod):
    if(mod == 'scRNASeq'):
        methods = evaluation_methods_dict['scRNASeq']
    elif(mod == 'CyTOF'):
        methods = evaluation_methods_dict['CyTOF']

    return(expand(output + 'rare-subtype-benchmarking/{{modality}}-Active-Learning-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{set}-{method}-predictions.tsv',
           set = list(range(training_set_AL_simulation)), method = methods))

def get_predictions(mod):
    if mod == 'scRNASeq':
        f = cell_type_predictions['scRNASeq']
    elif mod == 'snRNASeq':
        f = cell_type_predictions['snRNASeq']
    elif mod == 'CyTOF':
        f = cell_type_predictions['CyTOF']
    
    return (f)

def gt_pred_pdf_output(mod, selection_expansion_dict = selection_expansion_dict, evaluation_methods_dict = evaluation_methods_dict, train_test_seeds = train_test_seeds, cell_numbers = cell_numbers):
    rand_0 = []
    pred_rand_0 = [expand('{modality}/sel_{selection_procedure}/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf',
                        modality = mod,
                        initial = selection_expansion_dict[select]['initial'],
                        selection_procedure = [select], 
                        strat = selection_expansion_dict[select]['strategy'],
                        AL_alg = selection_expansion_dict[select]['AL_alg'],
                        rand = 0,
                        corruption = selection_expansion_dict[select]['corruption'],
                        knn = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        method = evaluation_methods_dict[mod],
                        s = train_test_seeds,
                        cell_num = cell_numbers) 
                        for select in selection_expansion_dict.keys()]
    
    for i in pred_rand_0:
        rand_0.extend(i)
    
    corr_0 = []
    pred_corr_0 = [expand('{modality}/sel_{selection_procedure}/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf',
                        modality = mod,
                        initial = selection_expansion_dict[select]['initial'],
                        selection_procedure = [select], 
                        strat = selection_expansion_dict[select]['strategy'],
                        AL_alg = selection_expansion_dict[select]['AL_alg'],
                        rand = selection_expansion_dict[select]['random_selection'],
                        corruption = 0,
                        knn = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        method = evaluation_methods_dict[mod],
                        s = train_test_seeds,
                        cell_num = cell_numbers) 
                        for select in selection_expansion_dict.keys()]

    for i in pred_corr_0:
        corr_0.extend(i)

    pred_lab = []
    pred_lab_temp = [expand('{modality}/sel_{selection_procedure}/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_lab_alg}-selection-{pred_lab_selection}-rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf',
                        modality = mod,
                        initial = selection_expansion_dict[select]['initial'],
                        selection_procedure = [select], 
                        strat = selection_expansion_dict[select]['strategy'],
                        AL_alg = selection_expansion_dict[select]['AL_alg'],
                        pred_lab_alg = ['rf', 'multinom'],
                        pred_lab_selection = ['entropy', 'maxp', '99p'], 
                        rand = 0,
                        corruption = 0,
                        knn = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        method = evaluation_methods_dict[mod],
                        s = train_test_seeds,
                        cell_num = cell_numbers) 
                        for select in selection_expansion_dict.keys()]

    for i in pred_lab_temp:
        pred_lab.extend(i)

    full = rand_0 + corr_0 + pred_lab
    
    return full

viz = {
    #'benchmark': expand(output + 'reports/overall-{modality}-{AL_alg}-benchmarking-cells.html', modality = modalities, AL_alg = AL_methods),
    # 'cross_cohort': output + 'reports/cross-cohort-eval.html',
    # 'random_selected_cells': expand(output + 'figures/selected_cells/{modality}/{selection}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA.pdf',
    #     modality = modalities, selection = 'random'),
    # 'seurat_selected_cells': expand(output + 'figures/selected_cells/{modality}/{selection}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{knn}-resolution-{res}.pdf',
    #     modality = modalities, selection = 'Seurat-clustering', knn = Seurat_neighbors, res = Seurat_resolution),
    # 'al_selected_cells_cor_entr': expand(output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-0-corr-all-knn_neighbors-{knn}-resolution-{res}.pdf',
    #     modality = modalities, selection = 'Active-Learning_entropy', strat = ['0.75_quant_entropy', '0.95_quant_entropy', 'highest_entropy'], al = AL_methods, knn = 'NA', res = 'NA'),
    # 'al_selected_cells_cor_maxp': expand(output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-0-corr-all-knn_neighbors-{knn}-resolution-{res}.pdf',
    #     modality = modalities, selection = 'Active-Learning_maxp', strat = ['0.05_quant_maxp', '0.25_quant_maxp', 'lowest_maxp'], al = AL_methods, knn = 'NA', res = 'NA'),
    # 'al_selected_cells_rand_entr': expand(output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-all-corr-0-knn_neighbors-{knn}-resolution-{res}.pdf',
    #     modality = modalities, selection = 'Active-Learning_entropy', strat = ['0.75_quant_entropy', '0.95_quant_entropy', 'highest_entropy'], al = AL_methods, knn = 'NA', res = 'NA'),
    # 'al_selected_cells_rand_maxp': expand(output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-all-corr-0-knn_neighbors-{knn}-resolution-{res}.pdf',
    #     modality = modalities, selection = 'Active-Learning_maxp', strat = ['0.05_quant_maxp', '0.25_quant_maxp', 'lowest_maxp'], al = AL_methods, knn = 'NA', res = 'NA'),
    # 'viz_corruption': expand(output + 'reports/corruption/vizualize-corruption-{modality}-cells-Init_{initial}-seed-{s}.html', 
    #     modality = modalities, initial = initial_selections, s = train_test_seeds),
    'gt_predicted_alluv': gt_pred_pdf_output('CyTOF') + gt_pred_pdf_output('scRNASeq') + gt_pred_pdf_output('snRNASeq')
}

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

rule visualize_selected_cells:
    input:
        sel_cells = expand(output + 'data/{{modality}}/{{selection}}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{cell_num}_cells.tsv', cell_num = cell_numbers, s = train_test_seeds)
    output:
        barplot = output + 'figures/selected_cells/{modality}/{selection}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-selected-cells.R'

rule visualize_selected_cells_AL_corrupt:
    input:
        sel_cells = expand(output + 'data/{{modality}}/{{selection}}/AL-batches-subset/Init_{initial}-strat-{{strat}}-ALAlg-{{al}}-rand_sel-0-corr-{cor}-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{cell_num}_cells.tsv', 
            initial = initial_selections, cell_num = cell_numbers, s = train_test_seeds, cor = corruption_percentages)
    params:
        cor_or_rand = "cor"
    output:
        barplot = output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-0-corr-all-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-selected-cells.R'

rule visualize_selected_cells_AL_rand:
    input:
        sel_cells = expand(output + 'data/{{modality}}/{{selection}}/AL-batches-subset/Init_{initial}-strat-{{strat}}-ALAlg-{{al}}-rand_sel-{rand}-corr-0-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{cell_num}_cells.tsv', 
            initial = initial_selections, cell_num = cell_numbers, s = train_test_seeds, rand = random_percentages)
    params:
        cor_or_rand = "rand"
    output:
        barplot = output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-all-corr-0-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-selected-cells.R'

rule visualize_corruption:
    input:
        AL_entropy = expand(output + 'data/{{modality}}/Active-Learning_entropy/AL-batches-subset/Init_{{initial}}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-0-corr-{corrupt}-knn_neighbors-NA-resolution-NA-seed-{{s}}-{cell_num}_cells.tsv',
            AL_alg = AL_methods, strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], corrupt = selection_expansion_dict['Active-Learning_entropy']['corruption'], cell_num = cell_numbers),
        AL_maxp = expand(output + 'data/{{modality}}/Active-Learning_maxp/AL-batches-subset/Init_{{initial}}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-0-corr-{corrupt}-knn_neighbors-NA-resolution-NA-seed-{{s}}-{cell_num}_cells.tsv',
            AL_alg = AL_methods, strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], corrupt = selection_expansion_dict['Active-Learning_maxp']['corruption'], cell_num = cell_numbers),
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    params:
        output_dir = output + 'reports/corruption/',
        AL_files_entropy = output + 'data/{modality}/Active-Learning_entropy/AL-batches-subset',
        AL_files_maxp = output + 'data/{modality}/Active-Learning_maxp/AL-batches-subset',
        pattern = 'Init_{initial}-*-rand_sel-0-*-seed-{s}-*.tsv'
    output:
        html = output + 'reports/corruption/vizualize-corruption-{modality}-cells-Init_{initial}-seed-{s}.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/ltri/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/visualize/visualize-corruption.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sce = '{input.sce}', AL_entropy_path = '{params.AL_files_entropy}', AL_maxp_path = '{params.AL_files_maxp}', modality = '{wildcards.modality}', pattern = '{params.pattern}'))\" "

def get_gt_predicted_input_files(mod, pred_labeling, imbalance = False):
    if pred_labeling:
        f = expand(output + 'predictive-labeling-benchmarking/{{modality}}/Init_{{initial}}-sel-{{selection_procedure}}-strat-{{strat}}-ALAlg-{{AL_alg}}-pred_alg-{{pred_lab_alg}}-rand_sel-{{rand}}-corr-{{corruption}}-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{method}-predictions-{cell_num}-cells-{{pred_lab_selection}}.tsv', 
                method = evaluation_methods_dict[mod],
                s = train_test_seeds,
                cell_num = cell_numbers)
    elif imbalance:
        f = expand(output + 'imbalance/rare-subtype-benchmarking/{{bal}}-{{similarity}}/Init-{{initial}}-{{modality}}-sel-{{selection_procedure}}-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corruption}}-knn_neighbors-{{knn}}-resolution-{{res}}-{method}-predictions-seed-{s}-{cell_num}-cells.tsv', 
                method = evaluation_methods_dict[mod],
                s = train_test_seeds,
                cell_num = 100)
    else:
        f = expand(output + 'rare-subtype-benchmarking/Init_{{initial}}-{{modality}}-sel_{{selection_procedure}}-strat-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corruption}}-knn_neighbors-{{knn}}-resolution-{{res}}-{method}-predictions-seed-{s}-{cell_num}-cells.tsv', 
                method = evaluation_methods_dict[mod],
                s = train_test_seeds,
                cell_num = cell_numbers)
    return f

rule compare_gt_to_predicted:
    input:
        sce = 'data/{modality}/{modality}-full.rds',
        predicted = lambda wildcards: get_gt_predicted_input_files(wildcards.modality, False)
    output:
        pdf = output + 'figures/gt_predictions/{modality}/sel_{selection_procedure}/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-gt-against-predicted.R'

rule compare_gt_to_predicted_pred_lab:
    input:
        sce = 'data/{modality}/{modality}-full.rds',
        predicted = lambda wildcards: get_gt_predicted_input_files(wildcards.modality, True)
    output:
        pdf = output + 'figures/gt_predictions_pred_labs/{modality}/sel_{selection_procedure}/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_lab_alg}-selection-{pred_lab_selection}-rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-gt-against-predicted.R'

# rule overall_benchmark:
#     input:
#         sce = 'data/{modality}/{modality}-full.rds',
#         predictions = lambda wildcards: get_predictions(wildcards.modality)
#     resources:
#         mem_mb=50000
#     log:
#         output + 'logs/benchmark-predictive-labeling-{modality}.log'
#     output:
#         acc = output + 'results/overall-{modality}-benchmarking-accuracies.tsv'
#     script:
#         'benchmarking/save-acc-overall-benchmarking.R'

# rule visualize_overall_benchmark:
#     input:
#         acc = output + 'results/overall-{modality}-benchmarking-accuracies.tsv'
#     output:
#         html = output + 'reports/overall-{modality}-{AL_alg}-benchmarking-cells.html'
#     params:
#         output_dir = output + 'reports/'
#     shell:
#         "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
#         "rmarkdown::render('pipeline/benchmarking/benchmarking.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
#         "params = list(acc = '{input.acc}', al_alg = '{wildcards.AL_alg}'))\" "
