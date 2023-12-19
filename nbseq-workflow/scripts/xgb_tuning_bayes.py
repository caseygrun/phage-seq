import warnings
from common import *

import pickle
import numpy as np

# monkey patch for skopt
np.int = int


import pandas as pd
import xgboost as xgb
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

warnings.filterwarnings("ignore", category=UserWarning,
                        message=r".*UserWarning: `use_label_encoder` is deprecated in 1.7.0 *")

with snakemake_log(snakemake) as logfile:

    print(f"Learning classifier for antigen '{antigen}' using design matrix '{dataset}'")
    print(f"Optimizing hyperparameters with Bayesian search")

    train_idx, X_train, X_unknown, y = get_design_for_ag(designs[dataset], ag_matrix, ag=antigen, shuffle=False)

    scale_pos_weight=(sum(y == 0) / sum(y > 0))

    estimator = xgb.XGBClassifier(learning_rate=0.1, n_estimators=100,
          max_depth=5, min_child_weight=1,
          gamma=0,
          subsample=0.8, colsample_bytree=0.8,
          objective= 'binary:logistic', eval_metric='auc',
          nthread=1, scale_pos_weight=scale_pos_weight,
          seed=1337, use_label_encoder=False)

    print("Base estimator: ")
    print(estimator)
    print()


    print("Search space:")
    search_spaces = {
        'max_depth': list(range(3,10,1)),
        'min_child_weight':list(range(1,6,1)),
        'subsample':       Real(low=0.5, high=1), #[0.1*i for i in range(6,10)],
        'colsample_bytree':Real(low=0.6, high=1), #[0.1*i for i in range(6,10)],
        'gamma':           Real(low=1e-5, high=100,  prior='log-uniform'),  #[i/10.0 for i in range(0,20)]
        'reg_alpha':       Real(low=1e-5, high=100, prior='log-uniform') # [1e-5, 1e-2, 0.1, 1, 100]
    }
    search = BayesSearchCV(
        estimator,
        search_spaces=search_spaces,
        cv=StratifiedShuffleSplit(n_splits=32, test_size=0.3),
        scoring='roc_auc',
        # n_iter=100,
        n_iter=50,
        n_points=snakemake.threads//2,
        n_jobs=snakemake.threads,
        verbose=2
    )
    print(search)
    print("Searching...")

    search.fit(X_train, y)

    print("Best parameters:")
    print(search.best_params_)

    print("Best score:")
    print(search.best_score_)

    # (summary, params, spaces)

    with open(snakemake.output[0], 'wb') as file:
        pickle.dump({ 'dataset':dataset, 'antigen':antigen, 'estimator': search }, file)
