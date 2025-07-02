rule vardict_paired_mode:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        region=config['resources'][genome_version]['GENOME_BED'],
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
    params:
        extra="",
        af_th=0.01,
    threads: 10
    conda:
        "vardict_java"
    shell:
        """
            vardict-java \
            -h \
            -th {threads} \
            -G {input.reference} \
            -N '{wildcards.sample}-T|{wildcards.sample}-NC' \
            -b '{input.Tum}|{input.normal}' \
            -Q 1 \
            -c 1 \
            -S 2 \
            -E 3 \
            -g 4 \
            {input.region} \
            | awk 'NR!=1' \
            | testsomatic.R \
            | var2vcf_paired.pl -N '{wildcards.sample}-T|{wildcards.sample}-NC' -f {params.af_th} > {output.vcf}
        """
        

rule vardict_filter_somatic:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict_filter.vcf"
    threads: 1
    params:
        caller='vardict'
    script:
        "../../../scripts/vcf_filter_somtic.R"

rule vardict_filter_germline:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/vardict_germline.vcf"
    threads: 1
    params:
        caller='vardict'
    script:
        "../../../scripts/vcf_filter_germline.R"