import warnings
from common import *

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import KFold, StratifiedKFold, RepeatedStratifiedKFold, GridSearchCV, StratifiedShuffleSplit

from nbseq.design import get_design_for_ag
from nbseq.predict import SerialHyperparamSearchResults, serial_hyperparam_search

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
    print(f"Optimizing hyperparameters with manual search")
    train_idx, X_train, X_unknown, y = get_design_for_ag(designs[dataset], ag_matrix, ag=antigen, shuffle=False)

    # (summary, params, spaces)
    out = serial_hyperparam_search(
        xgb.XGBClassifier(),
        X_train, y,

        starting_params=dict(
            learning_rate=0.1, n_estimators=100, max_depth=5,
            min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8
            ),
        param_grids=[
            {
                'max_depth':range(3,10,1),
                'min_child_weight':range(1,6,1)
            },
            {
                'subsample':        np.linspace(0.6,1,12), #range(0.6, 1, 0.05), #[0.1*i for i in range(6,10)],
                'colsample_bytree': np.linspace(0.6,1,12), #range(0.6, 1, 0.05), #[0.1*i for i in range(6,10)]
            },
            {
                'gamma':[i/10.0 for i in range(0,20)]
            },
            {
                'reg_alpha':[1e-5, 1e-2, 0.1, 1, 100]
            }
        ],
        constant_params=dict(
            objective= 'binary:logistic', eval_metric='auc',
            scale_pos_weight=(np.sum(y == 0)/np.sum(y > 0)),
            nthread=1, seed=1337, use_label_encoder=False,
        ),
        cv=StratifiedShuffleSplit(n_splits=32, test_size=0.3),
        plot=False,
        plot_progress=False,
        n_jobs=snakemake.threads
    )

    with open(snakemake.output[0], 'wb') as file:
        pickle.dump({ 'dataset':dataset, 'antigen':antigen, **out._asdict() }, file)
