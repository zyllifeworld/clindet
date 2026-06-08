rule octopus_call_paired:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/octopus/octopus.vcf"
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','octopus','sif'],image_url = config['singularity']['octopus']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.octopus.benchmark.txt"
    shell:
        """
        octopus \
            -R {input.ref} \
            -I {input.Tum} {input.NC} \
            --normal-sample {wildcards.sample}_NC \
            -t {input.region} \
            -o {output.vcf} \
            --threads {threads}
        """


rule octopus_paired_somatic:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/octopus/octopus.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/octopus.vcf"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools view \
            -i 'INFO/SOMATIC=1' \
            {input.vcf} \
            -Ov -o {output.vcf}
        """

rule octopus_paired_germline:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/octopus/octopus.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/octopus.vcf"
    conda:
        flexible_conda_env(config,['conda','clindet_main'],env_yaml = 'envs/clindet.yaml')
    shell:
        """
        bcftools view \
            -i 'INFO/SOMATIC!=1 && ((FILTER="PASS" && INFO/DP>30) || (FILTER="." && INFO/DP>30))'  \
            {input.vcf} -Ou | \
            bcftools view \
            -i 'INFO/SOMATIC!=1 && GT[0]!="0/0" && GT[0]!="0|0"' \
            -Ov -o {output.vcf}
        """

  