rule unpaired_call_variants_HaplotypeCaller:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/HaplotypeCaller.vcf",
    params:
        temp_directory=config['params']['java']['temp_directory']
    threads:10
    singularity:
        flexible_container_img(config,['singularity','gatk4','sif'],image_url = config['singularity']['gatk4']['repo'])
    benchmark:
        "{project}/{genome_version}/results/benchmarks/mut/{sample}.HaplotypeCaller.unpaired.benchmark.txt"
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && gatk \
        HaplotypeCaller -R {input.ref} \
        -I {input.Tum} \
        -O {output.vcf} \
        --intervals {input.bed} \
        --native-pair-hmm-threads {threads} --annotate-with-num-discovered-alleles -A UniqueAltReadCount -A ReferenceBases \
        -A PossibleDeNovo -A Coverage -A DepthPerAlleleBySample -A DepthPerSampleHC -A StrandBiasBySample -A StrandOddsRatio
        """