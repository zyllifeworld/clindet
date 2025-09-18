
rule unpaired_mutect2_call:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed,
        pon=config['resources'][genome_version]['WES_PON'],
        germ_vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Mutect2.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
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
        --intervals {input.bed}
        """

rule M2_filter_unpaired:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed,
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Mutect2.vcf",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/Mutect2_filter.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        FilterMutectCalls \
        -R {input.ref} \
        -V {input.vcf} \
        --stats {input.vcf}.stats \
        -O {output.vcf}
        """
