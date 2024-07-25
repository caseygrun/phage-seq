from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO
import hashlib
import pathlib

def md5(s):
	return hashlib.md5(s.encode('utf-8')).hexdigest()


def ungapped(seq, gap='-'):
    return seq.replace(gap, '')

def len_ungapped(seq, gap='-'):
    return len(seq) - seq.count(gap)

def feature_table_to_unique_sequences(feature_table):
	seqs = feature_table.index.unique().values
	print("Unique sequences:")
	print(seqs)
	seq_records = (SeqRecord(Seq(seq), id = md5(seq)) for seq in seqs)
	return seq_records

def seqrecords_to_dataframe(records):
	sequences = []
	ids = []
	names = []
	descriptions = []
	for record in records:
		ids.append(record.id)
		sequences.append(str(record.seq))
		names.append(record.name)
		descriptions.append(record.description)
	return pd.DataFrame({
		'sequence': sequences,
		'id': ids,
		'name': names,
		'description':descriptions
	})

def dataframe_to_seqrecords(df, seq_col = None, name_col = None, description_col=None, **kwargs):
	if name_col is None: names = df.index
	else: names = df.loc[:,name_col]

	if seq_col is None: seqs = df.iloc[:,0]
	else: seqs = df.loc[:,seq_col]

	if description_col is None: descs = [''] * len(seqs)
	else: descs = df.loc[:,description_col]

	seq_records = (SeqRecord(Seq(seq), id = name, description=desc, **kwargs) for seq,name,desc in zip(seqs,names, descs))
	return seq_records

def dataframe_to_fasta(df, file, **kwargs):
	return Bio.SeqIO.write(dataframe_to_seqrecords(df, **kwargs), file, format='fasta' )

def run_command(cmd, output_fp, env=None, **kwargs):
	print(' '.join(cmd))
	with open(output_fp, 'w') as output_f:
		subprocess.run(cmd, stdout=output_f, check=True, env=env, **kwargs)

def mkdirp_file(path):
	p = pathlib.Path(path).resolve().parent
	p.mkdir(parents=True, exist_ok=True)

def mkdirp(path):
	p = pathlib.Path(path)
	p.mkdir(parents=True, exist_ok=True)
