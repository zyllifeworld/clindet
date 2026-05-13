rule CNA_exomedepth:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        bed=get_sample_bed,
        ref=config['resources'][genome_version]['REFFA']
    output:
        rds="{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}/{sample}_exomedepth.rds",
        tsv="{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}/{sample}_exomedepth.tsv"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/exomedepth/{sample}",
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    benchmark:
        "{project}/{genome_version}/results/benchmarks/cnv/{sample}.exomedepth.benchmark.txt"
    conda: 
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    script:
        "../../../../scripts/ExomeDepth.R"