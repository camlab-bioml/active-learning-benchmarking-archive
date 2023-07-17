import numpy as np
import pandas as pd
import random

from sklearn.svm import SVC

import joblib


# ## Data processing
expression = pd.read_csv(snakemake.input['train'], sep = "\t", header = 0)

# read in markers and subset to those
with open(snakemake.input['markers'], 'r') as file:
  cell_types = yaml.safe_load(file)['cell_types']

m = list([list(x.values()) for x in cell_types.values()])
m = list(chain(*[x for x in m]))
m = list(chain(*[x for x in m]))
m = m + ['cell_type', 'cell_id']

expression = expression[m]


annotation = pd.read_csv(str(snakemake.input['annotation']), sep = '\t', header = 0)

if 'cell_selection' in dict(snakemake.wildcards).keys():
    expression = pd.merge(expression, annotation, on = 'cell_id')
    X_train = expression.drop(['cell_type_y', 'entropy', 'labeling', 'cell_type_x', 'cell_id', 'corrupted_cell_type', 'iteration', 'method', 'cell_num', 'params'], axis = 1)
    y_train = expression['cell_type_y']
else:
    expression = pd.merge(expression, annotation, on = 'cell_id')

    if snakemake.wildcards['selection_procedure'] == 'Active-Learning_entropy' or snakemake.wildcards['selection_procedure'] == 'Active-Learning_maxp':
        expression = expression.drop(['iteration', 'cell_num', 'corrupted_cell_type'], axis = 1)

    if snakemake.wildcards['selection_procedure'] == 'random':
        expression = expression.drop(['set','params'], axis = 1)

    if snakemake.wildcards['selection_procedure'] == 'NoMarkerSeurat-clustering' or snakemake.wildcards['selection_procedure'] == 'MarkerSeurat-clustering':
        expression = expression.drop(['params'], axis = 1)

    X_train = expression.drop(['cell_type_y', 'method', 'cell_type_x', 'cell_id'], axis = 1)
    y_train = expression['cell_type_y']

# ## ML Pipeline start
RSEED = 42

Classifier = SVC(probability=True, kernel='linear', random_state=0)
Classifier.fit(X_train, y_train)

joblib.dump(Classifier, snakemake.output['model'], compress = 1)

print(f"Training score: {grid.score(X_train, y_train)}")
print(f"The best parameters are: {grid.best_params_}")