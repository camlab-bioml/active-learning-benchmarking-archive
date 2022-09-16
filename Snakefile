
import pandas as pd
import numpy as np
import yaml
from itertools import chain

configfile: 'config/config.yml'
output = 'output/' + config['version'] + '/'

AL_methods = ['multinom', 'rf']
selection_procedures = ['random', 'Seurat-clustering', 'Active-Learning_entropy', 'Active-Learning_maxp']
evaluation_methods_dict = {
    'scRNASeq': ['scmap-cluster', 'scmap-sc', 'singleR', 'Random-Forest'], 
    'CyTOF': ['CyTOF-LDA', 'Random-Forest'],
    'snRNASeq': ['singleR', 'Random-Forest']
    }
evaluation_methods = list(chain.from_iterable(list(evaluation_methods_dict.values())))
train_test_seeds = list(range(10))
annotators = ['GroundTruth']
modalities = ['scRNASeq', 'CyTOF', 'snRNASeq']
data_splits = ['train', 'test']

# Seurat params
Seurat_neighbors = [10,20,30]
Seurat_resolution = [0.4,0.8,1.2]

random_sets = ['set1']
cell_numbers = [100, 250, 500, 900]
corruption_percentages = [0, 0.1, 0.2, 0.3, 1]


# Get markers for each cohort
with open(r'markers/scRNASeq.yml') as file:
    cell_types = yaml.full_load(file)

with open(r'markers/CyTOF.yml') as file:
    cytof_types = yaml.full_load(file)

with open(r'markers/snRNASeq.yml') as file:
    sn_cell_types = yaml.full_load(file)

# clean up marker names
scRNA_cell_types = list(cell_types['cell_types'].keys())
scRNA_cell_types_clean = []
for i in scRNA_cell_types:
    j = i.replace(' ', '-')
    scRNA_cell_types_clean.append(j)

cytof_cell_types = list(cytof_types['cell_types'].keys())
cytof_cell_types_clean = []
for i in cytof_cell_types:
    j = i.replace(' ', '-')
    cytof_cell_types_clean.append(j)

snRNASeq_cell_types = list(sn_cell_types['cell_types'].keys())
snRNASeq_cell_types_clean = []
for i in snRNASeq_cell_types:
    j = i.replace(' ', '-')
    snRNASeq_cell_types_clean.append(j)

all_cell_types = {
    'scRNASeq': scRNA_cell_types_clean, 
    'CyTOF': cytof_cell_types_clean, 
    'snRNASeq': snRNASeq_cell_types_clean
}

# Needed for cell type removal in AL - cell type ID's need to match exactly
original_cell_types = {
    'scRNASeq': scRNA_cell_types, 
    'CyTOF': cytof_cell_types, 
    'snRNASeq': snRNASeq_cell_types
}

include: 'pipeline/process-data.smk'
include: 'pipeline/cell-type-predictions.smk'
include: 'pipeline/simulate-active-learning.smk'
include: 'pipeline/visualizations.smk'
#include: 'pipeline/predictive-labeling.smk'
include: 'pipeline/imbalance.smk'
include: 'pipeline/rem-cell-type.smk'

rule all:
    input:
        process_data_output.values(),
        #cell_type_predictions.values(),
        #active_learner.values(),
        #viz.values(),
        #pred_lab.values(),
        #imbalance.values(),
        #rem_cell_type.values()

