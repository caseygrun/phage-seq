from common import snakemake_log
with snakemake_log(snakemake) as logfile:

    import nbseq.ft

    # hacky monkey patch because deicode refers to the deprecated np.int
    # https://numpy.org/devdocs/release/1.20.0-notes.html#using-the-aliases-of-builtin-types-like-np-int-is-deprecated
    import numpy as np
    np.int = int

    from deicode.rpca import auto_rpca

    ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])

    ordination, distance = auto_rpca(ft)

    distance.write(snakemake.output['distance_matrix'])
    ordination.write(snakemake.output['pcoa'])
