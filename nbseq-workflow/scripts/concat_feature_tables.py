import nbseq.ft

fts = [nbseq.ft.read_feature_table_biom(f) for f in snakemake.input]
ft = fts[0].concat(fts[1:])
nbseq.ft.write_feature_table_biom(ft, snakemake.output[0])
