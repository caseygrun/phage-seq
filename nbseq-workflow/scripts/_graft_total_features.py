

import nbseq.ft
import pandas as pd
import argparse

# python workflow/scripts/_graft_total_features.py --base results/tables/summary_total_features.txt \
# --summaries intermediate/{20220111,20211229}/{aa,cdrs,cdr3,align,aa_cluster,aa_filter}/summary_total_features.txt \
# --fts aa=intermediate/20211229/aa/feature_table.biom aa=intermediate/20220111/aa/feature_table.biom \
# cdrs=intermediate/20211229/cdrs/feature_table.biom cdrs=intermediate/20220111/cdrs/feature_table.biom \
# cdr3=intermediate/20211229/cdr3/feature_table.biom cdr3=intermediate/20220111/cdr3/feature_table.biom \
# align=results/runs/20211229/align/feature_table.biom align=results/runs/20220111/align/feature_table.biom \
# aa_cluster=intermediate/20211229/aa_cluster/feature_table.biom aa_cluster=intermediate/20220111/aa_cluster/feature_table.biom \
# aa_filter=intermediate/20211229/aa_filter/feature_table.biom aa_filter=intermediate/20220111/aa_filter/feature_table.biom \
# --output intermediate/summary_total_features.txt

parser = argparse.ArgumentParser(description='Replace columns of a summary file.')
parser.add_argument('--base', help='summary file to start with')
parser.add_argument('--row', default=None, 
                    help='Change row names to this')
parser.add_argument('--fts', dest='fts', nargs='+', action='extend',
                    help='column=path pairs for each feature table you want to use')
parser.add_argument('--summaries', dest='summaries', nargs='+', action='extend',
                    help='path to each summary file you want to replace')
parser.add_argument('--output', help='where to write the summary file')

args = parser.parse_args()

dfs = {}
output = pd.read_csv(args.base, index_col=0, sep="\t")
output.index = output.index.astype(str)

print(f"Starting from {args.base}:")
print(output)
print()

# PER-RUN VALUES
row_name = args.row
print(f"Add per-run summaries from {args.summaries}:")

for f in args.summaries:

    print(f"- Reading from {f}:")

    df = pd.read_csv(f, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if row_name is not None:
        df.rename(columns={df.columns[0]: row_name}, inplace=True)

    print(df)

    colname = df.columns[0]
    if colname in dfs:
        dfs[colname] = pd.concat([dfs[colname], df])
    else:
        dfs[colname] = df

print(f"Updating columns: {list(dfs.keys())}")
for colname, df in dfs.items():
    for idx in df.index.values:
        output.loc[idx, colname] = df.loc[idx, colname]

print(output)
print()
print()
del dfs

# TOTALS
print("Updating totals...")
print(f"Add TOTALS:")
column_fts = {}
for pair in args.fts:
    column, f = pair.split('=')
    if not (column in column_fts):
        column_fts[column] = []
    column_fts[column].append(f)

print(column_fts)
print()

for column, fts in column_fts.items():
    print(f"Updating column `{column}` from {fts}:")
    ids = set()
    for path in fts:
        print(f"- Reading from {path}:")
        ft = nbseq.ft.read_feature_table_biom(path)
        ids.update(ids, ft.ids('observation'))
        print(len(ids))
        del ft
    output.loc['TOTAL',column] = len(ids)
    del ids

print(f"Writing output to {args.output}:")
print(output)



output.to_csv(args.output, sep="\t", index=True)
