include: 'workflow/utils/functions.smk'
include: 'workflow/setup/setup.smk'

if config['run_type'] == 'wes':
    include: 'wrapper/wes.smk'
    rule all:
        input:
            rules.wes_pipeline.input
elif config['run_type'] == 'wgs':
    include: 'wrapper/wgs.smk'
    rule all:
        input:
            rules.wgs_pipeline.input
elif config['run_type'] == 'rna':
    include: 'wrapper/rna.smk'
    rule all:
        input:
            rules.rna_pipeline.input
elif config['run_type'] == 'build_b37':
    rule all:
        input:
            rules.build_b37_ref.input
elif config['run_type'] == 'build_hg38':
    rule all:
        input:
            rules.build_b37_ref.input
elif config['run_type'] == 'build_hg38':
    rule all:
        input:
            rules.build_b37_ref.input