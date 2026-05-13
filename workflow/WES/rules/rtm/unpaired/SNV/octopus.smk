rule octopus_call_tumor_only:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/octopus/octopus.vcf"
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','octopus','sif'],image_url = config['singularity']['octopus']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.octopus.benchmark.txt"
    shell:
        """

        """