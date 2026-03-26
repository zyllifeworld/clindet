# I can't test AmpliconSuite Singularity image on HPC,because AssertionError: Singularity version  3.5.3-1.1.el7 is not supported. Please upgrade to version 3.6 or higher.
# use conda env
def ampliconsuite_ref_version(wildcards):
    if wildcards.sample == 'b37':
        version = '--ref  GRCh37'
    elif wildcards.sample == 'hg19':
        version = '--ref  hg19'
    elif wildcards.sample == 'GRCh38':
        version = '--ref  GRCh38'
    else:
        version = ''
    return version

rule ampliconsuite:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        "{project}/{genome_version}/results/cnv/paired/ecDNA/{sample}"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/ecDNA/{sample}",
        ref_version=ampliconsuite_ref_version
    threads: 20
    conda: config['conda']['ampliconsuite']
    shell:
        """
        AmpliconSuite-pipeline.py 
        -o {params.wd} -s {wildcards.sample} \
        -t {threads} \
        --bam {input.Tum} \
        --normal_bam {input.NC} \
        {params.ref_version} --run_AA --run_AC
        """