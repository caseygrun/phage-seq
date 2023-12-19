import numpy as np
import nbseq.ft

ft = nbseq.ft.read_feature_table_biom(snakemake.input[0])

transform_functions = {
	'log1p': np.log1p,
	'sqrt': np.sqrt
}
transform_f = transform_functions[snakemake.wildcards['transform']]

ft = nbseq.ft.transform(ft, f=transform_f)

nbseq.ft.write_feature_table_biom(ft, snakemake.output[0])
