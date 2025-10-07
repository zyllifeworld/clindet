rule mutect2_call:
    input:
        Tum="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['resources'][genome_version]['REFFA'],
        pon=config['resources'][genome_version]['WES_PON'],
        germ_vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        bed=config['resources'][genome_version]['WES_BED'],
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
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
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2.vcf",
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/Mutect2_filter.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
    threads: 1
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        FilterMutectCalls \
        -R {input.ref} \
        -V {input.vcf} \
        --stats {input.vcf}.stats \
        -O {output.vcf}
        """