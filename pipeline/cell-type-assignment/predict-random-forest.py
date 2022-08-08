import numpy as np
import pandas as pd 

import joblib


## Load everything 
model = joblib.load(open(snakemake.input['model'], 'rb'))
expression = pd.read_csv(snakemake.input['test'], header = 0, sep = '\t')

true_labels = expression['cell_type']
cell_ids = expression['cell_id']
expression = expression.drop(['cell_type', 'cell_id'], axis = 1)

predicted = model.predict(expression)

output = pd.DataFrame({'cell_id': cell_ids,
                       'predicted_cell_type': predicted,
                       'prediction_params': 'Random-Forest-iterations_set-' + snakemake.wildcards['set'] + '-knn-' + snakemake.wildcards['neighbors'] + '-res-' + snakemake.wildcards['res'] + '-cell_numbers-' + snakemake.wildcards['cell_num'] + '-randomSelection-' + snakemake.wildcards['rand'] + '-corrupted-' + snakemake.wildcards['corrupt'],
                       'selection_procedure': snakemake.wildcards['selection_procedure'] + '-strategy-' + snakemake.wildcards['strat'],
                       'training_annotator': snakemake.wildcards['annotator'],
                       'modality': snakemake.params['modality']})

if 'cell_selection' in dict(snakemake.wildcards).keys():
    output['cell_selection'] = snakemake.wildcards['cell_selection']
else:
    output['cell_selection'] = 'NA'

output.to_csv(snakemake.output['predictions'], sep = '\t', index = False)