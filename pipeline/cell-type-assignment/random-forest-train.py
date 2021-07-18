import numpy as np
import pandas as pd
import random

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.decomposition import PCA
import joblib


# ## Data processing
expression = pd.read_csv(snakemake.input['train'], sep = "\t", header = 0)
X_train = expression.drop(['cell_type', 'cell_id'], axis = 1)
y_train = expression['cell_type']

# ## ML Pipeline start
RSEED = 42

steps = [('Scale', StandardScaler()), ('pca', PCA()), ('RandomForest', RandomForestClassifier())]
pipeline = Pipeline(steps)

parameters = {
    'pca__n_components': [100, 450],
    'RandomForest__max_features': [2, 4, 6, 10],
    'RandomForest__n_estimators': [50, 100, 150]
}

print("starting grid search")
grid = GridSearchCV(pipeline, param_grid=parameters, cv=10, return_train_score=True)

grid.fit(X_train, y_train)

joblib.dump(grid.best_estimator_, snakemake.output['model'], compress = 1)

print(f"Training score: {grid.score(X_train, y_train)}")
print(f"The best parameters are: {grid.best_params_}")