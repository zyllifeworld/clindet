rule unpaired_call_config_strelka:
    input:
        Tum="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=config['resources'][genome_version]['WES_BED'],
    output:
        dir=directory("{project}/{genome_version}/results/mut/vcf/{sample}/Manta"),
        tamp="{project}/{genome_version}/results/mut/vcf/{sample}/{sample}-Manta.log",
        indelcd="{project}/{genome_version}/results/mut/vcf/{sample}/Manta/results/variants/candidateSmallIndels.vcf.gz"
    params:
    conda:
        "strelka"
    shell:
        """
            [ ! -f {input.bed}.gz ] && gzip -k {input.bed}
            configManta.py \
            --tumorBam {input.Tum} \
            --runDir {output.dir} \
            --exome \
            --callRegions {input.bed}.gz \
            --referenceFasta {input.ref}
            {output.dir}/runWorkflow.py
            touch {output.tamp}
        """

rule unpaired_call_strelka_manta:
    input:
        Tum="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
        ref=config['resources'][genome_version]['REFFA'],
        indelcd="{project}/{genome_version}/results/mut/vcf/{sample}/Manta/results/variants/candidateSmallIndels.vcf.gz",
        bed=config['resources'][genome_version]['WES_BED'],
    output:
        dir=directory("{project}/{genome_version}/results/mut/vcf/{sample}/Strelka"),
        tamp="{project}/{genome_version}/results/mut/vcf/{sample}/{sample}-Strelka.log"
    params:
    conda:
        "strelka"
    shell:
        """
        configureStrelkaGermlineWorkflow.py \
        --bam {input.Tum} \
        --referenceFasta {input.ref} \
        --indelCandidates {input.indelcd} \
        --exome \
        --runDir {output.dir} \
        --callRegions {input.bed}.gz
        {output.dir}/runWorkflow.py -m local 
        touch {output.tamp}
        """



rule unpaired_strelka_filter:
    input:
        tamp="{project}/{genome_version}/results/mut/vcf/{sample}/{sample}-Strelka.log"
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/strelka.vcf"
    params:
        snp="{project}/{genome_version}/results/mut/vcf/{sample}/Strelka/results/variants/variants.vcf.gz"
    shell:
        """
        bcftools view -f 'PASS' {params.snp} -Ov -o {output.vcf}
        """

