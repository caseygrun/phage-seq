import Bio.SeqIO
import Bio.Align
import Bio.Align.substitution_matrices

def _initialize_alignment(library, reference_path, reference_length=101, reference_frame_start=0, peptide_length_threshold=50):
	"""
	perform a peptide alignment between each translated read and the reference peptide sequence, in order to identify CDRs
	"""
	global aligner
	global reference_translated_seq

	reference = next(Bio.SeqIO.parse(reference_path,"fasta"))
	reference_translated = reference[reference_frame_start:(reference_frame_start+reference_length*3)].translate(to_stop=True)
	reference_translated_seq = reference_translated.seq

	aligner = Bio.Align.PairwiseAligner()

	# ALPACA
	if library == 'alpaca':
		aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')

		# don't penalize matches to X, otherwise CDRs align badly
		aligner.substitution_matrix[:,'X'] = 0.5

		# the internal regions of the target (reference) sequence should be
		# generously long enough to accomodate variable length query sequences.
		aligner.target_internal_gap_score = -10000

		# in particular, there should
		aligner.target_left_gap_score = -10000

		# some reads may be longer than the target sequence; that's OK
		aligner.target_right_gap_score = 0

		# still need to penalize gaps somewhat, else
		aligner.query_open_gap_score = -1
		aligner.query_extend_gap_score = -0.5
		aligner.mismatch_score = -1

	# SYNTHETIC
	elif library == 'synthetic':

		# these seem to work better with default settings
		aligner.target_open_gap_score = -10000
		aligner.query_left_gap_score = -1000
		aligner.query_right_gap_score = -10



def _align_reference(seq):
	alignments = aligner.align(reference_translated_seq,seq)
	if len(alignments) < 1000:
		alignments = sorted(alignments)
#             [print(a) for a in sorted(alignments)]
	alignment = alignments[0]
	s = str(alignment).split("\n")[2]
	return s


def pairwise_align_reference_parallel(seqs, processes=None, **kwargs):
	import multiprocessing

	# # table = alignment_table[alignment_table[f'{library}_aligned'] & (alignment_table[f'{library}_length'] > peptide_length_threshold)]
	#
	# return table[f'{library}_translation'].apply(align_reference).rename(f'{library}_translation_align')

	with multiprocessing.Pool(processes, initializer = _initialize_alignment, initargs = kwargs) as pool:
		# TODO: consider using imap_unordered and then sorting at the end?
		aligned = pool.map(seqs, _align_reference)
	return dict(zip(seqs, aligned))

def pairwise_align_reference(seqs, **kwargs):
	_initialize_alignment(**kwargs)
	aligned = map(seqs, _align_reference)
	return dict(zip(seqs, aligned))
