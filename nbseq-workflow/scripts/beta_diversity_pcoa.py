from common import snakemake_log
with snakemake_log(snakemake) as logfile:

    import skbio
    import skbio.stats.ordination
    
    print(f"Reading input from {snakemake.input['distance_matrix']}")
    dm = skbio.DistanceMatrix.read(snakemake.input['distance_matrix'])

    print("Performing PCOA")
    pcoa = skbio.stats.ordination.pcoa(dm)

    print(f"Writing output to {snakemake.output['pcoa']}")
    pcoa.write(snakemake.output['pcoa'])
