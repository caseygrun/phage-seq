import re
import locale
import os

import pandas as pd

locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' )


pattern = re.compile(r"Pairs written \(passing filters\):\s+([\d,]+)\s+\(([\d\.]+)%\)")
files = []
reads = []
percent = []
for file in snakemake.input:
    with open(file) as f:
        for i, line in enumerate(f):
            match = re.match(pattern, line)
            if match is not None:
                rs, ps = match.groups()
                reads.append(locale.atoi(rs))
                percent.append(float(ps))
                files.append(os.path.splitext(os.path.basename(file))[0])
                break

out = pd.DataFrame({'sample':files,'reads_trimmed':reads,'percent_trimmed':percent})
out.to_csv(snakemake.output[0], index=False, sep="\t")
