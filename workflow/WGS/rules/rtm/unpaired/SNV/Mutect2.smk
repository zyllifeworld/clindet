rule unpaired_mutect2_call:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        germ_vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Mutect2_raw.vcf"
    params:
        temp_directory=config['params']['java']['temp_directory'],
        pon=lambda wildcards: get_config_value(config,['resources', wildcards.genome_version, 'WES_PON'], params='-pon', default=""),
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.mutect2.benchmark.txt"
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        Mutect2 -R {input.ref} \
        --native-pair-hmm-threads {threads} \
        -I {input.Tum} \
        -O {output.vcf} \
        {params.pon} \
        --germline-resource {input.germ_vcf}
        """

rule M2_filter_unpaired:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Mutect2_raw.vcf",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Mutect2.vcf"
    params:
        temp_directory=config['params']['java']['temp_directory'],
    threads: 10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        FilterMutectCalls \
        -R {input.ref} \
        -V {input.vcf} \
        --stats {input.vcf}.stats \
        -O {output.vcf}
        """
