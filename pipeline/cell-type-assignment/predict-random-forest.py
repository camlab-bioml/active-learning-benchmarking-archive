import numpy as np
import pandas as pd 

import joblib


## Load everything 
model = joblib.load(open(snakemake.input['model'], 'rb'))
expression = pd.read_csv(snakemake.input['test'])

true_labels = expression['cell_type']
cell_ids = expression['cell_id']
expression = expression.drop(['cell_type', 'cell_id'])

predicted = model.predict(expression)

output = pd.DataFrame({'cell_id': cell_ids,
                       'predicted_cell_type': predicted,
                       'prediction_params': NA,
                       'selection_procedure': snakemake.wildcards['selection_procedure'],
                       'training_annotator': snakemake.wildcards['training_annotator'],
                       'modality': snakemake.wildcards['modality']})

pd.to_csv(output, snakemake.output['predictions'])