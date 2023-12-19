import subprocess
import tempfile
from pathlib import Path
import pandas as pd
import nbseq.ft
from Bio import SeqIO

from common import snakemake_log
# def parse_id(id, seq, name, desc):
# 	id, *_ = id.split(";size=")
# 	return (id, seq, name, desc)

with open(snakemake.log[0],'w') as logfile: #snakemake_log(snakemake) as logfile:

	ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])
	pairs = pd.read_csv(snakemake.input['unique_pairs'], index_col=0)

	print("Constructing merged pairs...", file=logfile)

	seqs = (pairs['fwd'] +
		'NNNNNNNNNN' +
		pairs['rev'].apply(nbseq.utils.rc))

	abundance = (nbseq.ft.summarize(ft,
			axis='observation',
			id_name='pair_id',
			summary_name='abundance'))
	abundance = (abundance
		.set_index('pair_id')['abundance']
		.reindex_like(seqs))

	df = pd.DataFrame({
		'seq': seqs,
		'abundance': abundance,
	}, index=seqs.index)
	df = nbseq.utils.index_to_column(df, name='name')

	df['header'] = df['name'] + ';size=' + df['abundance'].astype(str)

	# with tempfile.TemporaryDirectory() as tmp:
	# 	tmp_pairs_fasta = Path(tmp) / 'input.fasta'
	# 	tmp_chimeras = Path(tmp) / 'chimeras.fasta'
	# 	tmp_nonchimeras = Path(tmp) / 'nonchimeras.fasta'
	# 	nbseq.utils.dataframe_to_fasta(df, tmp_pairs_fasta,
	# 		name_col='header')
	#
	# 	print("Running `vsearch`: ", file=logfile)
	# 	cmd = ['vsearch','--uchime_denovo',
	# 		str(tmp_pairs_fasta),
	# 		'--chimeras', str(tmp_chimeras),
	# 		'--nonchimeras', str(tmp_nonchimeras),
	# 		'--xsize' # don't print abundance information with chimeras
	# 	]
	# 	print(" ".join(cmd)+"", file=logfile)
	# 	subprocess.call(cmd, stdout=logfile, stderr=logfile)
	#
	# 	print("Done. Reading results...", file=logfile)
	# 	nonchimeras = nbseq.utils.seqrecords_to_dataframe(SeqIO.parse(tmp_nonchimeras,'fasta'))
	# 	chimeras = nbseq.utils.seqrecords_to_dataframe(SeqIO.parse(chimeras,'fasta'))
	# chimeras = pairs.loc[chimeras.ID,:]
	# nonchimeras = pairs.loc[nonchimeras.ID,:]
	nonchimeras = pairs
	chimeras = pd.DataFrame(data=None, columns=pairs.columns)


	print("Writing output...", file=logfile)

	# write chimera and non-chimera sequences back to csv files
	# if 'chimeras' in snakemake.output:
	chimeras.to_csv(snakemake.output['chimeras'])
	nonchimeras.to_csv(snakemake.output['unique_pairs'])

	# write filtered feature table
	# print("Writing filtered feature table...", file=logfile)
	ft = ft.filter(nonchimeras.index.values,axis='observation')
	nbseq.ft.write_feature_table_biom(ft, snakemake.output['feature_table'])
