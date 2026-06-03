include: 'workflow/utils/functions.smk'
include: 'workflow/setup/setup.smk'

import os

JAVA_TEMP_DIR = config.get("params", {}).get("java", {}).get("temp_directory", "tmp/java")
os.makedirs(JAVA_TEMP_DIR, exist_ok=True)

VALID_RUN_TYPES = [
    "wes",
    "wgs",
    "rna",
    "build_b37",
    "build_hg38",
]

RUN_TYPE = config.get("run_type", None)

if RUN_TYPE is None or RUN_TYPE == "":
    raise ValueError(
        "Missing required config parameter: run_type.\n"
        f"Please set config['run_type'] to one of: {', '.join(VALID_RUN_TYPES)}\n"
        "Example in config.yaml:\n"
        "run_type: wes"
    )

if RUN_TYPE not in VALID_RUN_TYPES:
    raise ValueError(
        f"Invalid config['run_type']: {RUN_TYPE}\n"
        f"Valid choices are: {', '.join(VALID_RUN_TYPES)}"
    )


if RUN_TYPE == 'wes':
    include: 'wrapper/wes.smk'
    rule all:
        input:
            rules.wes_pipeline.input

elif RUN_TYPE == 'wgs':
    include: 'wrapper/wgs.smk'
    rule all:
        input:
            rules.wgs_pipeline.input

elif RUN_TYPE == 'rna':
    include: 'wrapper/rna.smk'
    rule all:
        input:
            rules.rna_pipeline.input

elif RUN_TYPE == 'build_b37':
    rule all:
        input:
            rules.build_b37_ref.input

elif RUN_TYPE == 'build_hg38':
    rule all:
        input:
            rules.build_hg38_ref.input