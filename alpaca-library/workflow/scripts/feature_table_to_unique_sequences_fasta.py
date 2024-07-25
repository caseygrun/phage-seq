import nbseq
import pandas as pd
import Bio.SeqIO

feature_table = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
seq_records = nbseq.feature_table_to_unique_sequences(feature_table)
Bio.SeqIO.write(seq_records, snakemake.output[0], "fasta")
