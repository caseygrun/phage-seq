import nbseq.ft

ft = nbseq.ft.read_feature_table_sparse(snakemake.input[0])
nbseq.ft.write_feature_table_biom(ft, snakemake.output[0])
