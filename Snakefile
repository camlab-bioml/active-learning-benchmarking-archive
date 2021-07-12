
import pandas as pd
import numpy as np

configfile: 'config/config.yml'
output = 'output/' + config['version'] + '/'

selection_procedures = ['random']
annotators = ['test']
modalities = ['scRNASeq']
data_splits = ['train', 'test']

include: 'pipeline/process-data.smk'
include: 'pipeline/cell-type-predictions.smk'


rule all:
    input:
        process_data_output.values(),
        cell_type_predictions.values()