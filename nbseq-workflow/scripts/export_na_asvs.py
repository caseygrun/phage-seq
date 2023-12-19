import pandas as pd
from nbseq.utils import len_ungapped, dataframe_to_fasta

dfs = []
for library_csv_path in snakemake.input:
    df = pd.read_csv(library_csv_path, index_col=0)
    dfs.append(df)

df = pd.concat(dfs)

# if an ASV aligned to multiple references, keep the longer one
# TODO: this might not be the best choice
df['len'] = df['seq'].apply(len_ungapped)
df.index.rename('ASVID',inplace=True)
df = df.reset_index()

duplicated_ASVs = df.loc[df.duplicated('ASVID'),:]


df = df.sort_values(['ASVID','len'])
df = df.drop_duplicates('ASVID', keep='first')

# write output
df.set_index('ASVID')['seq'].to_csv(snakemake.output['csv'])
dataframe_to_fasta(df, snakemake.output['fasta'], seq_col='seq', name_col='ASVID')
