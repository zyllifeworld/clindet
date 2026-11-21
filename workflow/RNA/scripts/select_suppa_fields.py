#!/usr/bin/env python3
"""Convert Salmon quant.sf into SUPPA2 isoform TPM format.

Salmon Name values from the GENCODE transcript FASTA look like:
    ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|...
SUPPA2 expects the transcript ID alone in column 1.
"""

import csv
from pathlib import Path

quant_path = Path(snakemake.input.quant)
out_path = Path(snakemake.output.tpm)

out_path.parent.mkdir(parents=True, exist_ok=True)

with quant_path.open(newline="") as in_handle, out_path.open("w", newline="") as out_handle:
    reader = csv.reader(in_handle, delimiter="\t")
    header = next(reader)
    name_idx = header.index("Name")
    tpm_idx = header.index("TPM")

    # SUPPA2 expects a header line with sample names only (one column shorter
    # than the data rows; there is no label for the transcript-ID column).
    out_handle.write(f"{snakemake.wildcards.sample}\n")

    for row in reader:
        transcript_id = row[name_idx].split("|", 1)[0]
        out_handle.write(f"{transcript_id}\t{row[tpm_idx]}\n")
