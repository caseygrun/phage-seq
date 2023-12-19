from common import *

import pickle
import numpy as np
import pandas as pd

from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MaxAbsScaler
from sklearn.model_selection import KFold, StratifiedKFold, RepeatedStratifiedKFold, StratifiedShuffleSplit, GridSearchCV

from nbseq.design import get_design_for_ag

from skopt import BayesSearchCV
from skopt.space import Real, Categorical, Integer

with open(snakemake.input['designs'], 'rb') as file:
    designs = pickle.load( file )

with open(snakemake.input['ag_matrix'], 'rb') as file:
    ag_matrix = pickle.load( file )

dataset = snakemake.wildcards['design']
antigen = snakemake.wildcards['antigen']


with snakemake_log(snakemake) as logfile:
    import sys
    print(sys.stdout)
    print(sys.stderr)

    print(f"Learning classifier for antigen '{antigen}' using design matrix '{dataset}'")
    print(f"Optimizing hyperparameters with Bayesian search")

    train_idx, X_train, X_unknown, y = get_design_for_ag(designs[dataset], ag_matrix, ag=antigen, shuffle=False)


    estimator = make_pipeline(
        MaxAbsScaler(), 
        LogisticRegression(
            penalty='elasticnet',
            # penalty='l1',
            max_iter=4000,
            tol=1e-3,
            solver='saga',
            # solver='liblinear',
            class_weight='balanced',
            random_state=1337
        )
    )

    print("Base estimator: ")
    print(estimator)
    print()


    print("Search space:")
    search_spaces = {
        'logisticregression__l1_ratio': Real(low=0, high=1),
        'logisticregression__C': Real(low=1e-5, high=1000, prior='log-uniform')
    }
    search = BayesSearchCV(
        estimator,
        search_spaces=search_spaces,
        cv=StratifiedShuffleSplit(n_splits=32, test_size=0.3),
        # n_iter=50,
        n_iter=100,
        scoring='roc_auc',
        n_points=snakemake.threads//2,
        n_jobs=snakemake.threads,
        verbose=3
    )
    print(search)
    print("Searching...")

    search.fit(X_train, y)

    i = search.best_index_

    print("Best params:")
    print(search.best_params_)
    print("Best Score:")
    print(f"{search.cv_results_['mean_test_score'][i]} ± {search.cv_results_['std_test_score'][i]}")

    # (summary, params, spaces)

    with open(snakemake.output[0], 'wb') as file:
        pickle.dump({ 'dataset':dataset, 'antigen':antigen, 'estimator': search }, file)
