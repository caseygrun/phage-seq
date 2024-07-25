import os
import pandas as pd
import nbseq

if len(snakemake.log) > 0:
    log = open(snakemake.log, w)
    sys.stdout = log
    sys.stderr = log

# read dataframe, index should be ASV name/sequence, column should be
# translated AA sequence (we ignore the name here)
df = pd.read_csv(snakemake.input['translation'], index_col=0)
seqs = df.iloc[:, 0].values

print(f"Sequences to align ({len(seqs)}):")
print(seqs)

# aligned is a dict of AA_sequence:aligned_AA_sequence
aligned = nbseq.align_translated_peptides(seqs,
    processes = snakemake.threads,
    library = snakemake.wildcards['library'],
    reference_path = os.path.abspath(snakemake.input['reference']),
    reference_length = snakemake.params.library['reference_length_aa'],
    reference_frame_start=snakemake.params.library['reference_frame_start_nt'],
    min_length = snakemake.params.library['min_aa_length'])

# construct new DataFrame, index is ASV name/sequence, column 'aligned' is the
# aligned AA sequence. drop any rows which didn't align or were too short (as
# indicated by a blank `aligned`)
df_out = pd.DataFrame({ 'aligned': (aligned[x] for x in seqs) },
    index = df.index)
df_out.query('aligned != ""').to_csv(snakemake.output[0])
