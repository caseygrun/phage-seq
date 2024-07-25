from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO
import hashlib
import pandas as pd

from .utils import *


try:
	import Bio.Data.CodonTable

	Bio.Data.CodonTable.register_ncbi_table(
		name="Amber suppressing",
		alt_name=None,
		id=999,
		table={
			"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
			"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
			"TAT": "Y", "TAC": "Y",             "TAG": "Q",   # noqa: E241
			"TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
			"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
			"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
			"CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
			"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
			"ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
			"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
			"AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
			"AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
			"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
			"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
			"GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
			"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
		},
		stop_codons=["TAA", "TGA"],
		start_codons=["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"]
	)
except ImportError:
	print("Warning: BioPython is not installed. Some nbseq functions may not work.")


def translate_sam_alignment(samfile=None, samfile_path=None, reference=None, reference_path=None, reference_frame_start=0, suppress_amber=True):
	"""Load alignment records from a samfile/bamfile, find the first codon of each read, and translate that read.

	Parameters
	----------
	samfile : pysam.AlignmentFile
	samfile_path : str
	reference : Bio.Seq
	reference_path : str
	reference_frame_start : int, default=0

	Load data from `samfile_path` for `library` into a DataFrame. Use reference sequence
	at `reference_path` (whose relevant reading frame starts at `reference_frame_start`)
	to find the first complete codon and translate each aligned sequence.
	"""
	import pysam


	if samfile is None:
		if samfile_path is None:
			raise Exception("Must specify either `samfile` or `samfile_path`")
		samfile = pysam.AlignmentFile(samfile_path, "rb")

	# read and translate reference sequence
	# if reference is None:
	# 	if reference_path is None:
	# 		raise Exception("Must specify either `reference` or `reference_path`")
	# 	reference = next(Bio.SeqIO.parse(reference_path,"fasta"))
	# reference_translated = reference[reference_frame_start:].translate(to_stop=True)
	reference_translated = translate_reference(reference = reference, reference_path = reference_path, reference_frame_start = reference_frame_start)

	# allocate empty DataFrame, one row per unique ASV
	alignment_table = pd.Series(index=[row.query_sequence for row in samfile.fetch()], name='translation')

	print(f"Translating {len(alignment_table)} sequences from samfile...")

	for row in samfile.fetch():

		# identify the first complete codon
		# the correct reading frame for the reference starts with base `reference_frame_start`

		# codon:       0     | 1     | 2
		# refer: 0 1 | 2 3 4 | 5 6 7 | 8
		#            | G C G | G C C | G C T
		# query: - -   - - G   G C C   G C T
		# reference_start = 4
		# reference_frame_start = 2
		# first_codon = (reference_start - reference_frame_start) // 3 + 1 = (4 - 2) // 3 + 1
		first_codon = (((row.reference_start - reference_frame_start) // 3) + 1)

		# codon:           | 0     | 1     | 2
		# refer: 0 1 2 3 4 | 5 6 7 | 8 9 10
		#                  | G C G | G C C | G C T
		# reference_start:       ^
		#
		# query:       0 1   2 3 4   5 6 7   8 9 10
		#        - - - A A   A A G   G C C   G C T
		# query_alignment_start: ^
		# want first_base =          ^
		#
		# reference_start = 7
		# reference_frame_start = 5
		# query_alignment_start = 4
		# first_codon = (reference_start - reference_frame_start) // 3 + 1 = (7 - 5) // 3 + 1 = 1
		# first_base = query_alignment_start + (first_codon * 3 - (reference_start - reference_frame_start)) = 4 + (1 * 3 - (7-5)) = 4 + 3 - 2 = 5
		first_query_base = row.query_alignment_start + (first_codon * 3 - (row.reference_start - reference_frame_start))

		# extract an in-frame sequence
		seq_from_first_full_codon = row.query_sequence[first_query_base:]
		seq = Seq(seq_from_first_full_codon)

		# translate
		# TODO: handle amber suppression
		if suppress_amber:
			translated = seq.translate(to_stop=True, table=999)
		else:
			translated = seq.translate(to_stop=True)

		# pad the translated sequence to be same length as the reference
		# add dashes to beginning of sequence if it starts within the reference
		translated_padded = ('-' * first_codon) + str(translated)
		# trim to no longer than the length of the reference sequence
		translated_padded = translated_padded[:len(reference_translated)]
		# add dashes to make same length as the reference sequence (if shorter)
		translated_padded = translated_padded + ('-' * (len(reference_translated) - len(translated_padded)))
		assert len(translated_padded) == len(reference_translated)

		# alignment_table.loc[row.seq,f'{library}_aligned'] = (not row.is_unmapped)
		# alignment_table.loc[row.seq,f'{library}_length']  = len(translated)
		#
		# if library == 'alpaca':
		# 	alignment_table.loc[row.seq,f'{library}_inframe'] = len(translated) > 100
		# else:
		# 	alignment_table.loc[row.seq,f'{library}_inframe'] = len(translated) > 70

		alignment_table.loc[row.seq] = translated_padded
		# alignment_table.loc[row.seq,f'{library}_padded'] = translated_padded
	return alignment_table


def translate_reference(reference = None, reference_path = None, reference_frame_start = 0):
	# read and translate reference sequence
	if reference is None:
		if reference_path is None:
			raise Exception("Must specify either `reference` or `reference_path`")
		reference = next(Bio.SeqIO.parse(reference_path,"fasta"))
	reference_translated = reference[reference_frame_start:].translate(to_stop=True)
	return reference_translated

# populate_alignment_table(
# 	'alpaca',
# 	Path(SCRATCH_PATH) / "output/align/ASVs/table-alpaca-sorted.bam",
# 	Path(SCRATCH_PATH) / "output/align/references/2019-02-14_Alpaca_VHH/2019-02-14_CG023_Alpaca_VHH_sequences_good_trimmed_consensus.fasta")
#
# populate_alignment_table(
# 	'synthetic',
# 	Path(SCRATCH_PATH) / "output/align/ASVs/table-synthetic-CDR3_18aa-sorted.bam",
# 	Path(SCRATCH_PATH) / "output/align/references/2020-03-26_Feimei-Liu_pLibF_nanobody_display_variable_CDR3/2020-03-26_Feimei-Liu_pLibF_nanobody_display_CDR3_18aa.fasta",
# 	reference_frame_start = 953)
#
# table_align = pd.merge(left=table_tidy, left_on='ASV', right=alignment_table, right_index=True)



# from Bio import Align
# from Bio.Align import substitution_matrices
#
#
# def peptide_alignment(alignment_table, library, reference_path, reference_length=101, reference_frame_start=0, peptide_length_threshold=50):
# 	"""
# 	perform a peptide alignment between each translated read and the reference peptide sequence, in order to identify CDRs
# 	"""
# 	reference = next(Bio.SeqIO.parse(reference_path,"fasta"))
# 	reference_translated = reference[reference_frame_start:(reference_frame_start+reference_length*3)].translate(to_stop=True)
# 	print(reference_translated.seq)
#
# 	aligner = Align.PairwiseAligner()
#
# 	# ALPACA
# 	if library == 'alpaca':
# 		aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
#
# 		# don't penalize matches to X, otherwise CDRs align badly
# 		aligner.substitution_matrix[:,'X'] = 0.5
#
# 		# the internal regions of the target (reference) sequence should be
# 		# generously long enough to accomodate variable length query sequences.
# 		aligner.target_internal_gap_score = -10000
#
# 		# in particular, there should
# 		aligner.target_left_gap_score = -10000
#
# 		# some reads may be longer than the target sequence; that's OK
# 		aligner.target_right_gap_score = 0
#
# 		# still need to penalize gaps somewhat, else
# 		aligner.query_open_gap_score = -1
# 		aligner.query_extend_gap_score = -0.5
# 		aligner.mismatch_score = -1
#
# 	# SYNTHETIC
# 	elif library == 'synthetic':
#
# 		# these seem to work better with default settings
# 		aligner.target_open_gap_score = -10000
# 		aligner.query_left_gap_score = -1000
# 		aligner.query_right_gap_score = -10
#
#
# 	print(aligner)
#
# 	def align_reference(seq):
# 		alignments = aligner.align(reference_translated.seq,seq)
# #         print(len(alignments))
# #         try:
# #         print(alignments[0])
# #         alignment = sorted(alignments)[0]
# 		if len(alignments) < 1000:
# 			alignments = sorted(alignments)
# #             [print(a) for a in sorted(alignments)]
# 		alignment = alignments[0]
# 		print(alignment)
# #         print("================================================================")
# #         except Exception as e:
# #             print(e)
# #             print("!!!! OVERFLOW")
# #             print()
# #             print(alignments[0])
# #             print()
# #             alignment = alignments[0]
#
# 		s = str(alignment).split("\n")[2]
# 		print()
# 		return s
#
# 	table = alignment_table[alignment_table[f'{library}_aligned'] & (alignment_table[f'{library}_length'] > peptide_length_threshold)]
#
# 	return table[f'{library}_translation'].apply(align_reference).rename(f'{library}_translation_align')


# alpaca_peptide_align = peptide_alignment(
# 	alignment_table,
# 	'alpaca',
# 	Path(SCRATCH_PATH) / "output/align/references/2019-02-14_Alpaca_VHH/2019-02-14_CG023_Alpaca_VHH_sequences_good_trimmed_consensus.fasta",
# 	reference_length = 135,
# 	# trim beginning of leader sequence to keep alignment from tripping over GG's
# 	reference_frame_start = 12*3,
# 	peptide_length_threshold = 100
# )
#
# synthetic_peptide_align = peptide_alignment(
# 	alignment_table,
# 	'synthetic',
# 	Path(SCRATCH_PATH) / "output/align/references/2020-03-26_Feimei-Liu_pLibF_nanobody_display_variable_CDR3/2020-03-26_Feimei-Liu_pLibF_nanobody_display_CDR3_18aa.fasta",
# 	reference_length = 100,
# 	reference_frame_start = 953
# )

def align_translated_peptides(seqs, **kwargs):
	from .peptide import pairwise_align_reference_parallel
	return pairwise_align_reference_parallel(seqs, **kwargs)

def trim_align(max_length):
	def do_trim(x):
		if pd.isna(x):
			return x
		else:
			x = x[:max_length]
			if len(x) < max_length: x = x + ('-'*(max_length - len(x)))
			return x
	return do_trim



CDR_positions = {
	'alpaca': {
		'CDR1': (17, 26),
		'CDR2': (46, 53),
		'CDR3': (94, 115)
	},
	'synthetic': {
		'CDR1': (0, 6),
		'CDR2': (23, 29),
		'CDR3': (70, 87),
	}
}

def extract_CDRs(seqs, CDRs=None, library=None, push_left=False):
	"""Extract sequences of CDRs at defined positions from a sequence of aligned peptide sequences

	Parameters
	----------
	seqs : pd.Series
		Translated and aligned peptide sequences
	CDRs : dict, optional
		keys are the names of the CDR(s), values are `(start, end)` position tuples
	library : str, optional
		if CDRs is omitted, this can be given instead and values will be taken from `CDR_positions`
	push_left : bool, Default=False
		if True, remove leading dashes ('-') from each CDR sequence, in order to push the aligned sequence to the left

	Returns
	-------
	pd.DataFrame
		one column for each key in `CDRs`, indexed the same as `seqs`

	"""
	def extract_subsequence(start,end):
		length = end - start
		def do_extract(seq):
			if pd.isna(seq): return seq
			else:
				if push_left:
					seq = seq[start:end].lstrip('-')
					return seq + ('-' * (length-len(seq)) )
				else: return seq[start:end]

		return do_extract

	if CDRs is None:
		CDRs = CDR_positions[library]

	CDR_df = pd.DataFrame(columns = CDRs.keys(), index = seqs.index)
	for subseq_name, start_end in CDRs.items():
		CDR_df[subseq_name] = seqs.apply(extract_subsequence(*start_end))

	return CDR_df


# alpaca_CDR_positions =
# table_align = extract_CDRs(table_align, 'alpaca', alpaca_CDR_positions)
#
# synthetic_CDR_positions =
#
# table_align = extract_CDRs(table_align, 'synthetic', synthetic_CDR_positions)


def segment_lens(table, segments, len_f=len_ungapped):
	"""Calculate the length of one or more segments

	table : pd.DataFrame
	segments : list of str
		list of column names; calculate the length of each of these columns with len_f
	len_f : callable
	"""
	return pd.concat([
		(table[segments]
                    .applymap(len_f, na_action='ignore')),
		table.drop(segments, axis='columns')
	], axis=1)
