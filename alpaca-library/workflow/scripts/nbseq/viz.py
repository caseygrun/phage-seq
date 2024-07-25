import weblogo
import numpy as np


# chemistry_ambig = weblogo.ColorScheme(
# 		[
# 			weblogo.SymbolColor("BSTYCZ", "green", "polar"),
# 			weblogo.SymbolColor("JNQ", "purple", "neutral"),
# 			weblogo.SymbolColor("KRH", "blue", "basic"),
# 			weblogo.SymbolColor("DE", "red", "acidic"),
# 			weblogo.SymbolColor("PAWFLIMV", "black", "hydrophobic"),
# 			weblogo.SymbolColor("-", "#DDDDDD", "gap")
# 		],
# 		alphabet=weblogo.protein_alphabet
# )
chemistry_ambig = weblogo.ColorScheme(
		[
			weblogo.SymbolColor("BSTNQYZ", "green", "polar"),
			weblogo.SymbolColor("JNGC", "purple", "neutral"),
			weblogo.SymbolColor("KRH", "blue", "basic"),
			weblogo.SymbolColor("DE", "red", "acidic"),
			weblogo.SymbolColor("AVILMFW", "black", "hydrophobic"),
			weblogo.SymbolColor("-", "#DDDDDD", "gap")
		],
		alphabet=weblogo.protein_alphabet
)

def seq_abundances_to_count_matrix(seqs, abundances, alphabet):
	if not alphabet:
		raise ValueError("No alphabet")

	if len(seqs) != len(abundances):
		raise ValueError("Must have exactly one count per sequence")

	# number of characters in alphabet
	N = len(alphabet)
	ords = [alphabet.ords(s) for s in seqs]

	# length of first (= each) sequence
	L = len(ords[0])

	#     counts = [[0, ] * N for l in range(0, L)]
	# pre-allocate
	counts = np.zeros((L,N))

	# for each sequence
	for o, count in zip(ords, abundances):

		# check length is correct
		if len(o) != L:
			raise ValueError("Sequences are of incommensurate lengths. Cannot tally.")

		# build up counts
		for j, n in enumerate(o):
			if n < N:
				counts[j,n] += count

	return counts

def seq_abundances_to_motif(seqs, abundances, alphabet):
	from weblogo.matrix import Motif

	counts = seq_abundances_to_count_matrix(seqs, abundances, alphabet)
	return Motif(alphabet, counts)

def make_seqlogo(
		table, seq_col='sequence', abundance_col='abundance',
		logo_start=1, logo_end=None, color_scheme=chemistry_ambig, **kwargs):

	alphabet = weblogo.protein_alphabet

	valid_rows = (~table[seq_col].isna() & (table[seq_col] != ''))
	sequences = table.loc[valid_rows, seq_col]

	if logo_end is None:
		logo_end = len(sequences[0])

	if abundance_col is not None:
		abundances = table.loc[valid_rows, 'abundance']
		counts = seq_abundances_to_motif(sequences, abundances, alphabet)
		logodata = weblogo.LogoData.from_counts(alphabet, counts)
	else:
		sequences = [weblogo.seq.Seq(row, name=ix, alphabet=alphabet)
			for (i, (ix, row)) in enumerate(sequences.iteritems())]
		logodata = weblogo.LogoData.from_seqs(
			weblogo.seq.SeqList(sequences, alphabet=alphabet))

	logooptions = weblogo.LogoOptions(
		color_scheme=color_scheme,
		logo_start=logo_start,
		logo_end=logo_end,
		stacks_per_line=logo_end,
		 **{'unit_name':'probability',
			'show_boxes':False,
			'show_errorbars':False,
			**kwargs})
	logoformat = weblogo.LogoFormat(logodata, logooptions)
	return (logodata, logoformat)

def display_logo(logo, format='png'):
	if format == 'png':
		from IPython.core.display import display_png
		display_png(weblogo.png_formatter(*logo),raw=True)
	elif format=='svg':
		from IPython.core.display import display_svg
		display_svg(weblogo.svg_formatter(*logo),raw=True)
	elif format=='pdf':
		from IPython.core.display import display_pdf
		display_pdf(weblogo.pdf_formatter(*logo),raw=True)
