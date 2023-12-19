import nbseq.select

from common import snakemake_log
with snakemake_log(snakemake) as logfile:
    nbseq.select.enrichment_lme4(
        df_path=snakemake.input[0],
        output_path=snakemake.output[0],
        ag_col=snakemake.wildcards['antigen'],
        threads=snakemake.threads)
