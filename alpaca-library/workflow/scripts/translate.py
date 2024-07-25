import pandas as pd
import nbseq


alignment_table = nbseq.translate_sam_alignment(samfile_path=snakemake.input['alignment'],
    reference_path=snakemake.input['reference'],
    reference_frame_start=snakemake.params.library['reference_frame_start_nt'])

alignment_table.to_frame().to_csv(snakemake.output[0])
