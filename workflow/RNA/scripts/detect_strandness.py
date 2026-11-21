#!/usr/bin/env python3
"""Resolve splicing strandness from Salmon's auto-detection.

Modes:
  auto                          read aux_info/meta_info.json
  fr-firststrand                use as-is
  fr-secondstrand               use as-is
  fr-unstranded                 use as-is
"""

import json
from pathlib import Path

quant_path = Path(snakemake.input.quant)
out_path = Path(snakemake.output.txt)
mode = snakemake.params.mode
default = snakemake.params.default

LIBRARY_TYPE_MAP = {
    "ISR": "fr-firststrand",
    "ISF": "fr-secondstrand",
    "IU": "fr-unstranded",
    "U": "fr-unstranded",
}

result = default
if mode and mode != "auto":
    result = mode
else:
    meta_path = quant_path.parent / "aux_info" / "meta_info.json"
    if meta_path.exists():
        data = json.loads(meta_path.read_text())
        for library_type in data.get("library_types") or []:
            if library_type in LIBRARY_TYPE_MAP:
                result = LIBRARY_TYPE_MAP[library_type]
                break

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(result + "\n")
