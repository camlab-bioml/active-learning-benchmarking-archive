
import pandas as pd
import numpy as np
import yaml

configfile: 'config/config.yml'
output = 'output/' + config['version'] + '/'

selection_procedures = ['random']
annotators = ['test']
modalities = ['scRNASeq', 'CyTOF']
data_splits = ['train', 'test']

with open(r'markers/scRNA.yml') as file:
    cell_types = yaml.full_load(file)

scRNA_cell_types = list(cell_types['cell_types'].keys())
scRNA_cell_types_clean = []
for i in scRNA_cell_types:
    j = i.replace(' ', '-')
    scRNA_cell_types_clean.append(j)

include: 'pipeline/process-data.smk'
include: 'pipeline/cell-type-predictions.smk'


rule all:
    input:
        process_data_output.values(),
        cell_type_predictions.values()