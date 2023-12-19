import pandas as pd
import argparse


# python workflow/scripts/_graft_summaries.py --base results/tables/sample_summary_reads.txt --summaries intermediate/{20220111,20211229}/{aa,cdrs,cdr3,align,aa_cluster,aa_filter}/summary_reads.txt --output intermediate/sample_summary_reads.txt
# python workflow/scripts/_graft_summaries.py --base results/tables/sample_summary_features.txt --summaries intermediate/{20220111,20211229}/{aa,cdrs,cdr3,align,aa_cluster,aa_filter}/summary_features.txt --output intermediate/sample_summary_features.txt

parser = argparse.ArgumentParser(description='Replace columns of a summary file.')
parser.add_argument('--base', help='summary file to start with')
parser.add_argument('--row', default=None, 
                    help='Change row names to this')
parser.add_argument('--summaries', dest='summaries', nargs='+', action='extend',
                    help='path to each summary file you want to replace')
parser.add_argument('--output', help='where to write the summary file')

args = parser.parse_args()

dfs = {}
output = pd.read_csv(args.base, index_col=0, sep="\t")
print(f"Starting from {args.base}:")
print(output)

row_name = args.row

print(f"Add columns from {args.summaries}:")

for f in args.summaries:

    print(f"- Reading from {f}:")

    df = pd.read_csv(f, sep="\t", index_col=0)
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
    output[colname] = df.reindex(output.index)

print(f"Writing output to {args.output}:")
print(output)

output.to_csv(args.output, sep="\t", index=True)

