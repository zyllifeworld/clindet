rule call_variants_HaplotypeCaller:
    input:
        # Tum="{project}/{genome_version}/results/recal/{sample}.bam",
        Tum="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/HaplotypeCaller.vcf",
    threads:10
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        bed=config['resources'][genome_version]['WES_BED']
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        HaplotypeCaller -R {input.ref} \
        --intervals {params.bed} \
        -I {input.Tum} \
        -O {output.vcf} \
        -DF MappingQualityAvailableReadFilter \
        --native-pair-hmm-threads {threads} --annotate-with-num-discovered-alleles -A UniqueAltReadCount -A ReferenceBases \
        -A PossibleDeNovo -A Coverage -A DepthPerAlleleBySample -A DepthPerSampleHC -A StrandBiasBySample -A StrandOddsRatio
        """

rule norm_filter_HaplotypeCaller:
    input:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/HaplotypeCaller.vcf",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/vcf_nf/{sample}/HaplotypeCaller.vcf",
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
    shell:
        """
        bcftools norm -f {input.ref} -m +both -Ov {input.vcf} | bcftools filter -e 'QUAL<20 | FORMAT/DP[0] < 10' > {output.vcf}
        """