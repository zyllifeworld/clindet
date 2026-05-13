rule mutect2_call:
    input:
        Tum="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['resources'][genome_version]['REFFA'],
        pon=config['resources'][genome_version]['WES_PON'],
        germ_vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2_raw.vcf"
    params:
        temp_directory=config['params']['java']['temp_directory'],
        bed=config['resources'][genome_version]['WES_BED'],
    threads: 10
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.mutect2_call.benchmark.txt"
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        Mutect2 -R {input.ref} \
        --native-pair-hmm-threads {threads} \
        -I {input.Tum} \
        -O {output.vcf} \
        -pon {input.pon} \
        --germline-resource {input.germ_vcf} \
        --intervals {params.bed}
        """

rule M2_filter_unpaired:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2_raw.vcf",
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2.vcf"
    params:
        temp_directory=config['params']['java']['temp_directory'],
    threads: 1
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