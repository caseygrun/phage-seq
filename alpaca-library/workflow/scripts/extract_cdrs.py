import nbseq
import pandas as pd

# read dataframe, index should be ASV name/sequence, column should be
# translated, aligned AA sequence (we ignore the name here)
df = pd.read_csv(snakemake.input[0], index_col=0)
seqs = df.iloc[:, 0]


cdr_positions = None
if snakemake.params.CDRs:
    print("CDR positions:")
    cdr_positions = snakemake.params.CDRs
else:
    print("Using default CDR positions")

cdrs = nbseq.extract_CDRs(seqs, library=snakemake.wildcards['library'], CDRs=cdr_positions)
cdrs.to_csv(snakemake.output[0])
