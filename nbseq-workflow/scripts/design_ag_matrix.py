from common import *

import numpy as np
import pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler, MaxAbsScaler, PowerTransformer
from nbseq.predict import synthesize_input_controls
from nbseq.ft import summarize_filtering
import nbseq
import nbseq.design
import nbseq.select
import nbseq.ft
from nbseq.utils import *
from anndata import AnnData
import pickle
from importlib import reload
import warnings

warnings.filterwarnings("ignore", category=FutureWarning,
                        message=r".*FutureWarning: pandas.Int64Index is deprecated.*")
warnings.filterwarnings("ignore", category=FutureWarning,
                        message=r".*FutureWarning: pandas.core.index is deprecated.*")

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    # from pandas import MultiIndex, Int64Index
    import xgboost as xgb

# shut up with the deprecation warnings already
# %env PYTHONWARNINGS = ignore:: FutureWarning

with snakemake_log(snakemake) as logfile:

    # 'intermediate/learning/designs.pickle'
    design_path = snakemake.output['designs']

    # 'intermediate/learning/ag_matrix.pickle'
    ag_matrix_path = snakemake.output['ag_matrix']

    # load experiment
    ex = nbseq.Experiment.from_files(tree_cdr3=None, tree_aa=None, fd_cdr3=None)

    ft_cdr3_high = nbseq.ft.filter_abundance_prevalence(
        ex.fts.cdr3, 
        min_abundance=get_snakemake_param(snakemake,'min_abundance_cdr3', 10), 
        min_prevalence=get_snakemake_param(snakemake,'min_prevalence_cdr3', 5)
    )
    ft_cdr3_high_r = nbseq.ft.to_relative(ft_cdr3_high)

    print("Filtering CDR3 feature table:")
    summarize_filtering(ex.fts.cdr3, ft_cdr3_high)

    ft_aa_high = nbseq.ft.filter_abundance_prevalence(
        ex.fts.aa, 
        min_abundance=get_snakemake_param(snakemake, 'min_abundance_aa', 5),
        min_prevalence=get_snakemake_param(snakemake, 'min_prevalence_aa', 3),
    )
    ft_aa_high_r = nbseq.ft.to_relative(ft_aa_high)

    print("Filtering AA feature table:")
    summarize_filtering(ex.fts.aa, ft_aa_high)

    print("Synthesizing input controls...")

    print("- CDR3: ")
    ft_cdr3_high_w_controls = nbseq.ft.query(synthesize_input_controls(ft_cdr3_high_r, ex.ag_names),
                                            "io == 'i' & kind == '+'", axis='obs')

    print("- AA: ")
    ft_aa_high_w_controls = nbseq.ft.query(synthesize_input_controls(ft_aa_high_r, ex.ag_names),
                                        "io == 'i' & kind == '+'", axis='obs')

    print("Assembling antigen matrix: ")
    ag_matrix = ft_cdr3_high_w_controls.obs.drop_duplicates(
        subset='name').set_index('name')[ex.ag_names]

    print(ag_matrix)

    print("Saving antigen matrix to " + ag_matrix_path)
    # assert len(design) == len(ag_matrix)

    with open(ag_matrix_path, 'wb') as file:
        pickle.dump(ag_matrix, file)

    print("Assembling design matrices...")

    enrichment_kwargs = dict(fillna='R3i')
    ft_cdr3_high_w_controls.var_names.name = 'CDR3ID'
    enrichment_cdr3 = nbseq.select.calculate_enrichment_ad(
        ft_cdr3_high_w_controls, log=True, **enrichment_kwargs)  # , comparison=('R1i','R5i'))

    print(
        f"Calculating enrichment statistics for space 'CDR3' with {enrichment_kwargs}:")
    nbseq.select.summarize_enrichment(enrichment_cdr3)

    ft_aa_high_w_controls.var_names.name = 'feature'

    enrichment_kwargs = dict(fillna='R3i')
    enrichment_aa = nbseq.select.calculate_enrichment_ad(
        ft_aa_high_w_controls, log=True, feature_col='feature', **enrichment_kwargs)  # , comparison=('R1i','R5i'))

    print(
        f"Calculating enrichment statistics for space 'AA' with {enrichment_kwargs}:")
    nbseq.select.summarize_enrichment(enrichment_aa)

    design_enr = nbseq.design.enrichment_ft_to_design(
        ft_enrichment=enrichment_cdr3, fillna=None)
    # design_enr_dense = nbseq.design.enrichment_ft_to_design(
    #     ft_enrichment=enrichment_cdr3, fillna=0)

    design_last = nbseq.design.last_round(ft_cdr3_high_w_controls)
    # design_last_dense = nbseq.design.last_round(ft_cdr3_high_w_controls, fillna=0)

    design_last_enr = nbseq.design.concat_designs([design_enr, design_last])
    # design_last_enr_dense = nbseq.design.concat_designs(
    #     [design_enr_dense, design_last_dense])

    design_abd = nbseq.design.concat_rounds_colwise(ft_cdr3_high_w_controls)
    # design_abd_dense = nbseq.design.concat_rounds_colwise(
    #     ft_cdr3_high_w_controls, fillna=0)

    ft_cdr3_high_w_controls_log = nbseq.ft.transform(
        ft_cdr3_high_w_controls, np.log10)
    ft_cdr3_high_w_controls_log.X = sparse_drop_na(ft_cdr3_high_w_controls_log.X)
    fillna_log_abundance = np.floor(np.min(ft_cdr3_high_w_controls_log.X.data))
    design_abd_log = nbseq.design.concat_rounds_colwise(
        ft_cdr3_high_w_controls_log)

    design = nbseq.design.concat_designs([design_abd,       design_enr])
    design_log = nbseq.design.concat_designs([design_abd_log,   design_enr])
    # design_dense = nbseq.design.concat_designs(
    #     [design_abd_dense, design_enr_dense])

    design_last_log = nbseq.design.last_round(ft_cdr3_high_w_controls_log)
    design_last_log_dense = nbseq.design.last_round(
        ft_cdr3_high_w_controls_log, fillna=fillna_log_abundance)

    design_last_log_enr = nbseq.design.concat_designs(
        [design_last_log,       design_enr])
    # design_last_log_enr_dense = nbseq.design.concat_designs(
    #     [design_last_log_dense, design_enr_dense])

    ft_aa_high_w_controls_log = nbseq.ft.transform(ft_aa_high_w_controls, np.log10)
    ft_aa_high_w_controls_log.X = sparse_drop_na(ft_aa_high_w_controls_log.X)
    fillna_log_abundance = np.floor(np.min(ft_aa_high_w_controls_log.X.data))

    design_aa_abd = nbseq.design.concat_rounds_colwise(ft_aa_high_w_controls)
    design_aa_abd_log = nbseq.design.concat_rounds_colwise(
        ft_aa_high_w_controls_log)

    design_aa_last = nbseq.design.last_round(
        ft_aa_high_w_controls, identifier='feature')
    design_aa_enr = nbseq.design.enrichment_ft_to_design(
        ft_enrichment=enrichment_aa, fillna=None)

    design_aa_last_enr = nbseq.design.concat_designs(
        [design_aa_last,       design_aa_enr])

    design_aa = nbseq.design.concat_designs([design_aa_abd,     design_aa_enr])
    design_aa_log = nbseq.design.concat_designs([design_aa_abd_log, design_aa_enr])


    # there appears to be some pathology with the yeo-johnson transform when a variable has large offset but low variance.
    # this appears to be the case for the enrichment variables. This appears to fix
    # https://github.com/scikit-learn/scikit-learn/issues/14959
    trans_center_yeojohnson = make_pipeline(
        StandardScaler(with_std=False),
        PowerTransformer(standardize=True),
    )
    # trans = PowerTransformer() #(method='box-cox')


    def standardize_df(df, transformer):
        X_trans = transformer.fit_transform(df.values)
        return pd.DataFrame(X_trans, index=df.index.copy(), columns=df.columns.copy())


    def standardize_ad(ad, transformer):
        X_trans = transformer.fit_transform(ad.X)
        ad = ad.copy()
        ad.X = X_trans
        return ad


    def standardize(X, transformer):
        if isinstance(X, pd.DataFrame):
            return standardize_df(X, transformer)
        else:
            return standardize_ad(X, transformer)


    designs = {
        'CDR3:enr':               design_enr,
        # 'CDR3:enr_dense':         design_enr_dense,
        # 'CDR3:enr_std-ma':        standardize(design_enr, MaxAbsScaler()),

        'CDR3:R5+enr':            design_last_enr,
        # 'CDR3:R5+enr_dense':      design_last_enr_dense,

        'CDR3:log(R5)+enr':         design_last_log_enr,
        # 'CDR3:log(R5)+enr_dense':   design_last_log_enr_dense,
        # 'CDR3:log(R5)+enr_std-c':   standardize(design_last_log_enr_dense.to_df(), StandardScaler(with_std=True)),
        # 'CDR3:log(R5)+enr_std-mm':  standardize(design_last_log_enr_dense.to_df(), MinMaxScaler()),
        # 'CDR3:log(R5)+enr_std-ma':  standardize(design_last_log_enr, MaxAbsScaler()),

        'CDR3:R2345+enr':             design,
        # 'CDR3:R2345+enr_dense':       design_dense,
    }

    # from joblib import Parallel, delayed
    # parallel = Parallel(n_jobs=-1)

    # _designs = {
    #     'CDR3:enr_std-cyj':          lambda: standardize(design_enr_dense.to_df(), trans_center_yeojohnson),
    #     'CDR3:R5+enr_std-cyj':       lambda: standardize(design_last_enr_dense.to_df(), trans_center_yeojohnson),
    #     'CDR3:log(R5)+enr_std-cyj':  lambda: standardize(design_last_log_enr_dense.to_df(), trans_center_yeojohnson),
    #     'CDR3:R2345+enr_std-cyj':    lambda: standardize(design_dense.to_df(), trans_center_yeojohnson),
    # }

    # designs.update(
    #     dict(parallel(
    #             delayed(lambda k, v: (k, v()))(k, v) for k, v in _designs.items()
    #     ))
    # )

    designs.update({
        'CDR3:log(R2345)+enr':        design_log,
        # 'CDR3:log(R2345)+enr_std-mm': standardize(design_log.to_df().fillna(0), MinMaxScaler()),
        # 'CDR3:log(R2345)+enr_std-ma': standardize(design_log, MaxAbsScaler()),
    })

    designs.update({
        'aa:enr': design_aa_enr,
        'aa:R5+enr': design_aa_enr,

        'aa:R2345+enr':             design_aa,
        'aa:log(R2345)+enr':        design_aa_log,
        # 'aa:log(R2345)+enr_std-ma': standardize(design_aa_log, MaxAbsScaler())
    })

    designs = {key: designs[key] for key in sorted(designs.keys())}

    print(f"Available design matrices: {list(designs.keys())}")

    print("Saving design matrices to "+design_path)

    with open(design_path, 'wb') as file:
        pickle.dump(designs, file)
