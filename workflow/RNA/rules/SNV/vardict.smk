rule unpaired_vardict_single_mode:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=config['resources'][genome_version]['WES_BED'],
        bam="{project}/{genome_version}/results/mut/dedup/{sample}.split.bam",
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/VarDict/{sample}.vardict.vcf"
    params:
        extra="",
        bed_columns="-c 1 -S 2 -E 3 -g 4",  # Optional, default is -c 1 -S 2 -E 3 -g 4
        allele_frequency_threshold="0.01",  # Optional, default is 0.01
        post_scripts="teststrandbias.R | var2vcf_valid.pl -A -N "
    threads: 10
    conda: "clindet"
    shell:
        """
        vardict-java -G {input.reference} \
        {params.extra} \
        -th {threads} \
        {params.bed_columns} \
        -N '{wildcards.sample}_T' \
        -b {input.bam} \
        {input.regions} | {params.post_scripts} '{wildcards.sample}_T' -E > {output.vcf}
        """

rule unpaired_filter_vardict:
    input:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/VarDict/{sample}.vardict.vcf"
    output:
        vcf="{project}/{genome_version}/results/mut/vcf/{sample}/vardict.vcf"
    shell:
        """
        bcftools view -i 'FILTER="PASS" && INFO/DP>=20 && INFO/VD>=5'  {input.vcf} -Ov -o {output.vcf}
        """
