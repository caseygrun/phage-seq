import pandas as pd
import nbseq.ft

input = snakemake.input

res = pd.DataFrame([], columns=input.keys(), index=['TOTAL'])
for column, paths in input.items():
    ids = set()
    for path in paths[:2]:
        ft = nbseq.ft.read_feature_table_biom(path)
        ids.update(ids, ft.ids('observation'))
        del ft
    res.loc['TOTAL',column] = len(ids)

res.to_csv(snakemake.output[0], sep="\t",index=True)
