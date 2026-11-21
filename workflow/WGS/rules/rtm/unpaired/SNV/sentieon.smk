rule call_variants_sentieon_unpaired:
    input:
        tumor="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config["resources"][genome_version]["REFFA"],
        pon=config["resources"][genome_version]["WGS_PON"],
        dbsnp=config["resources"][genome_version]["DBSNP"]
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/sentieon.tnscope.pre_filter.vcf.gz"
    params:
        outdir="{project}/{genome_version}/results/vcf/unpaired/{sample}"
    threads: 10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.sentieon_tnscope_unpaired.benchmark.txt"
    shell:
        """
        {config[softwares][sentieon][call]} driver \
            -t {threads} \
            -r {input.ref} \
            -i {input.tumor} \
            --algo TNscope \
                --tumor_sample {wildcards.sample}_T \
                --dbsnp {input.dbsnp} \
                --pon {input.pon} \
                {output.vcf}
        """

rule filter_sentieon_tumor_only:
    input:
        vcf="{project}/{genome_version}/results/vcf/tumor_only/{sample}/sentieon.tnscope.pre_filter.vcf.gz"
    output:
        vcf="{project}/{genome_version}/results/vcf/tumor_only/{sample}/sentieon.vcf.gz",
        tbi="{project}/{genome_version}/results/vcf/tumor_only/{sample}/sentieon.vcf.gz.tbi"
    params:
        filter_script=config["softwares"]["sentieon"]["tnscope_filter"]
    threads: 2
    conda:
        flexible_conda_env(config,["conda", "clindet_main"],env_yaml="envs/clindet.yaml")
    shell:
        """
        {config[softwares][sentieon][call]} pyexec {params.filter_script} \
            -v {input.vcf} \
            --tumor_sample {wildcards.sample}_T \
            {output.vcf}

        tabix -p vcf {output.vcf}
        """
