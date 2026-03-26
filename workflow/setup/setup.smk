rule all:
    input:
        rules.build_reference.input

include: "rules/env_setup.smk"
include: "rules/human_b37.smk"