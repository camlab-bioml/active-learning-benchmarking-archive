
rem_cell_type_reports = []
rctreport = [expand(output + 'reports/rem_cell_type/Rem-{rem_cell_type}-{modality}-{{strat}}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}.html',
    rem_cell_type = original_cell_types[modality],
    modality = [modality],
    AL_alg = AL_methods,
    rand = 0,
    corrupt = 0)
    for modality in original_cell_types.keys()]

for element in rctreport:
    rem_cell_type_reports.extend(element)

rem_cell_type = {
    'reports': expand(rem_cell_type_reports, strat = ['highest_entropy', 'lowest_maxp'])
}

rule rem_cell_type:
    input:
        markers = 'markers/{modality}.yml',
        sce = 'data/{modality}/{modality}-train.rds'
    params:
        min_cells = 10,
        output_dir = output + 'reports/rem_cell_type/'
    output:
        html = output + 'reports/rem_cell_type/Rem-{rem_cell_type}-{modality}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/removed-cell-types/rem-cell-type.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(markers = '{input.markers}', sce = '{input.sce}', needed_cells = '{params.min_cells}', rem_cell_type = '{wildcards.rem_cell_type}', "
        "al_alg = '{wildcards.AL_alg}', strat = '{wildcards.strat}', rand = '{wildcards.rand}'))\" "